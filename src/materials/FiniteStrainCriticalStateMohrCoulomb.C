/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "FiniteStrainCriticalStateMohrCoulomb.h"
#include <math.h> // for M_PI

template<>
InputParameters validParams<FiniteStrainCriticalStateMohrCoulomb>()
{
  InputParameters params = validParams<FiniteStrainPlasticBase2>();
  params.addRequiredRangeCheckedParam<Real>("mc_cohesion", "mc_cohesion>=0", "Mohr-Coulomb cohesion");
  params.addRequiredRangeCheckedParam<Real>("mc_friction_angle", "mc_friction_angle>=0 & mc_friction_angle<=60", "Mohr-Coulomb friction angle in degrees");
  params.addParam<Real>("mc_cohesion_residual", "Mohr-Coulomb cohesion at infinite hardening.  If not given, this defaults to mc_cohesion, ie, perfect plasticity");
  params.addParam<Real>("mc_friction_angle_residual", "Mohr-Coulomb friction angle in degrees at infinite hardening.  If not given, this defaults to mc_friction_angle, ie, perfect plasticity");
  params.addRangeCheckedParam<Real>("mc_cohesion_rate", 0, "mc_cohesion_rate>=0", "Cohesion = mc_cohesion_residual + (mc_cohesion - mc_cohesion_residual)*exp(-mc_cohesion_rate*plasticstrain).  Set to zero for perfect plasticity");
  params.addRangeCheckedParam<Real>("mc_friction_angle_rate", 0, "mc_friction_angle_rate>=0", "friction_angle = mc_friction_angle_residual + (mc_friction_angle - mc_friction_angle_residual)*exp(-mc_friction_angle_rate*plasticstrain).  Set to zero for perfect plasticity");
  params.addRangeCheckedParam<Real>("mc_edge_smoother", 25.0, "mc_edge_smoother>=0 & mc_edge_smoother<=30", "Smoothing parameter: the edges of the cone are smoothed by the given amount.");
  params.addRequiredRangeCheckedParam<Real>("mc_tip_smoother", "mc_tip_smoother>=0", "Smoothing parameter: the cone vertex at mean = cohesion*cot(friction_angle), will be smoothed by the given amount.  Typical value is 0.1*cohesion");
  params.addParam<Real>("mc_lode_cutoff", "If the second invariant of stress is less than this amount, the Lode angle is assumed to be zero.  This is to gaurd against precision-loss problems, and this parameter should be set small.  Default = 0.00001*((yield_Function_tolerance)^2)");
  params.addRequiredParam<Real>("initial_void_ratio", "initial void ratio");
  params.addRequiredParam<Real>("gamma", "intercept of CSL in e-p' space");
  params.addRequiredParam<Real>("lambda", "slope of CSL in e-p' space");
  params.addRequiredParam<Real>("xi", "exponential of CSL in e-p' space");
  params.addRequiredParam<Real>("p_ref", "reference stress corresponding to 100 kPa");
  params.addRequiredParam<Real>("psi_0", "constant to adjust dilatancy angle");
  params.addRequiredParam<bool>("fixed_elastic_moduli", "If elastic moduli are not fixed, they will be computed using c_g and poisson's ratio");
  params.addParam<Real>("c_g", "constant for shear modulus");
  params.addParam<Real>("poisson_ratio", "poisson's ratio");
  params.addClassDescription("Non-associative Mohr-Coulomb plasticity with hardening/softening");

  return params;
}

FiniteStrainCriticalStateMohrCoulomb::FiniteStrainCriticalStateMohrCoulomb(const InputParameters & parameters) :
    FiniteStrainPlasticBase2(parameters),
    _initial_void_ratio(getParam<Real>("initial_void_ratio")),
    _gamma(getParam<Real>("gamma")),
    _lambda(getParam<Real>("lambda")),
    _xi(getParam<Real>("xi")),
    _p_ref(getParam<Real>("p_ref")),
    _psi_0(getParam<Real>("psi_0")),
    _fixed_elastic_moduli(getParam<bool>("fixed_elastic_moduli")),
    _c_g(getParam<Real>("c_g")),
    _poisson_ratio(getParam<Real>("poisson_ratio")),
    _cohesion(getParam<Real>("mc_cohesion")),
    _phi(getParam<Real>("mc_friction_angle")*M_PI/180.0),
    _cohesion_residual(parameters.isParamValid("mc_cohesion_residual") ? getParam<Real>("mc_cohesion_residual") : _cohesion),
    _phi_residual(parameters.isParamValid("mc_friction_angle_residual") ? getParam<Real>("mc_friction_angle_residual")*M_PI/180.0 : _phi),
    _cohesion_rate(getParam<Real>("mc_cohesion_rate")),
    _phi_rate(getParam<Real>("mc_friction_angle_rate")),
    _small_smoother2(std::pow(getParam<Real>("mc_tip_smoother"), 2)),
    _tt(getParam<Real>("mc_edge_smoother")*M_PI/180.0),
    _costt(std::cos(_tt)),
    _sintt(std::sin(_tt)),
    _cos3tt(std::cos(3*_tt)),
    _sin3tt(std::sin(3*_tt)),
    _cos6tt(std::cos(6*_tt)),
    _sin6tt(std::sin(6*_tt)),
    _lode_cutoff(parameters.isParamValid("mc_lode_cutoff") ? getParam<Real>("mc_lode_cutoff") : 1.0E-5*std::pow(_f_tol[0], 2)),
    _mc_internal(declareProperty<Real>("mc_internal")),
    _mc_internal_old(declarePropertyOld<Real>("mc_internal")),
    _mc_max_principal(declareProperty<Real>("mc_max_principal_stress")),
    _mc_min_principal(declareProperty<Real>("mc_min_principal_stress")),
    _yf(declareProperty<Real>("mc_yield_function")),
    _void_ratio(declareProperty<Real>("mc_void_ratio")),
    _dilatancy_angle(declareProperty<Real>("mc_dilatancy_angle"))
{
  if (_lode_cutoff < 0)
    mooseError("mc_lode_cutoff must not be negative");

  // check Abbo et al's convexity constraint (Eqn c.18 in their paper)
  Real sin_angle = std::sin(std::max(_phi, _phi_residual));
  Real rhs = std::sqrt(3)*(35*std::sin(_tt) + 14*std::sin(5*_tt) - 5*std::sin(7*_tt))/16/std::pow(std::cos(_tt), 5)/(11 - 10*std::cos(2*_tt));
  if (rhs <= sin_angle)
    mooseError("Mohr-Coulomb edge smoothing angle is too small and a non-convex yield surface will result.  Please choose a larger value");
}

void FiniteStrainCriticalStateMohrCoulomb::initQpStatefulProperties()
{
  _mc_internal[_qp] = 0;
  _mc_internal_old[_qp] = 0;
  _mc_max_principal[_qp] = 0;
  _mc_min_principal[_qp] = 0;
  _yf[_qp] = 0.0;
  _void_ratio[_qp] = _initial_void_ratio;
  FiniteStrainPlasticBase2::initQpStatefulProperties();
}

void
FiniteStrainCriticalStateMohrCoulomb::postReturnMap()
{
  // Record the value of the yield function
  _yf[_qp] = _f[_qp][0];

  // Record the value of the internal parameter
  _mc_internal[_qp] = _intnl[_qp][0];

  // update void ratio
  _void_ratio[_qp] = _initial_void_ratio + (_deformation_gradient[_qp].det() - 1) * (_initial_void_ratio + 1);

  // set min and max void ratio
  if (_void_ratio[_qp] < 0.6)
    _void_ratio[_qp] = 0.6;
  else if (_void_ratio[_qp] > 1.0)
    _void_ratio[_qp] = 1.0;

  // set dilatancy angle
  Real p = -1.0 * _stress[_qp].trace()/3.0;
  p = std::max(p, 0.0);
  _dilatancy_angle[_qp] = _psi_0 * (_gamma - _lambda*std::pow(p/_p_ref, _xi) - _void_ratio[_qp]) * 180.0 / M_PI;

  // update elasticity tensor if elastic moduli are not fixed
  if (_fixed_elastic_moduli == false)
    _elasticity_tensor[_qp] = get_elasticity_tensor(_stress[_qp], _void_ratio[_qp]);

  // Record the maximum principal stress
  std::vector<Real> eigvals;
  _stress[_qp].symmetricEigenvalues(eigvals);
  _mc_max_principal[_qp] = eigvals[2];
  _mc_min_principal[_qp] = eigvals[0];
}



unsigned int
FiniteStrainCriticalStateMohrCoulomb::numberOfInternalParameters()
{
  return 1;
}

void
FiniteStrainCriticalStateMohrCoulomb::yieldFunction(const RankTwoTensor &stress, const std::vector<Real> & intnl, std::vector<Real> & f)
{
  Real mean_stress = stress.trace()/3.0;
  Real sinphi = std::sin(phi(intnl[0]));
  Real cosphi = std::cos(phi(intnl[0]));
  Real sin3Lode = stress.sin3Lode(_lode_cutoff, 0);
  if (std::abs(sin3Lode) <= _sin3tt)
  {
    // the non-edge-smoothed version
    std::vector<Real> eigvals;
    stress.symmetricEigenvalues(eigvals);
    f.assign(1, mean_stress*sinphi + std::sqrt(_small_smoother2 + 0.25*std::pow(eigvals[2] - eigvals[0] + (eigvals[2] + eigvals[0] - 2*mean_stress)*sinphi, 2)) - cohesion(intnl[0])*cosphi);
  }
  else
  {
    // the edge-smoothed version
    Real aaa, bbb, ccc;
    abbo(sin3Lode, sinphi, aaa, bbb, ccc);
    Real kk = aaa + bbb*sin3Lode + ccc*std::pow(sin3Lode, 2);
    Real sibar2 = stress.secondInvariant();
    f.assign(1, mean_stress*sinphi + std::sqrt(_small_smoother2 + sibar2*std::pow(kk, 2)) - cohesion(intnl[0])*cosphi);
  }
}

RankTwoTensor
FiniteStrainCriticalStateMohrCoulomb::df_dsig(const RankTwoTensor & stress, const Real sin_angle)
{
  Real mean_stress = stress.trace()/3.0;
  RankTwoTensor dmean_stress = stress.dtrace()/3.0;
  Real sin3Lode = stress.sin3Lode(_lode_cutoff, 0);
  if (std::abs(sin3Lode) <= _sin3tt)
  {
    // the non-edge-smoothed version
    std::vector<Real> eigvals;
    std::vector<RankTwoTensor> deigvals;
    stress.dsymmetricEigenvalues(eigvals, deigvals);
    Real tmp = eigvals[2] - eigvals[0] + (eigvals[2] + eigvals[0] - 2*mean_stress)*sin_angle;
    RankTwoTensor dtmp = deigvals[2] - deigvals[0] + (deigvals[2] + deigvals[0] - 2*dmean_stress)*sin_angle;
    Real denom = std::sqrt(_small_smoother2 + 0.25*std::pow(tmp, 2));
    return dmean_stress*sin_angle + 0.25*tmp*dtmp/denom;
  }
  else
  {
    // the edge-smoothed version
    Real aaa, bbb, ccc;
    abbo(sin3Lode, sin_angle, aaa, bbb, ccc);
    Real kk = aaa + bbb*sin3Lode + ccc*std::pow(sin3Lode, 2);
    RankTwoTensor dkk = (bbb + 2*ccc*sin3Lode)*stress.dsin3Lode(_lode_cutoff);
    Real sibar2 = stress.secondInvariant();
    RankTwoTensor dsibar2 = stress.dsecondInvariant();
    Real denom = std::sqrt(_small_smoother2 + sibar2*std::pow(kk, 2));
    return dmean_stress*sin_angle + (0.5*dsibar2*std::pow(kk, 2) + sibar2*kk*dkk)/denom;
  }
}

void
FiniteStrainCriticalStateMohrCoulomb::dyieldFunction_dstress(const RankTwoTensor & stress, const std::vector<Real> & intnl, std::vector<RankTwoTensor> & df_dstress)
{
  Real sinphi = std::sin(phi(intnl[0]));
  df_dstress.assign(1, df_dsig(stress, sinphi));
}


void
FiniteStrainCriticalStateMohrCoulomb::dyieldFunction_dintnl(const RankTwoTensor & stress, const std::vector<Real> & intnl, std::vector<std::vector<Real> > & df_dintnl)
{
  df_dintnl.resize(1);

  Real sin_angle = std::sin(phi(intnl[0]));
  Real cos_angle = std::cos(phi(intnl[0]));
  Real dsin_angle = cos_angle*dphi(intnl[0]);
  Real dcos_angle = -sin_angle*dphi(intnl[0]);

  Real mean_stress = stress.trace()/3.0;
  Real sin3Lode = stress.sin3Lode(_lode_cutoff, 0);
  if (std::abs(sin3Lode) <= _sin3tt)
  {
    // the non-edge-smoothed version
    std::vector<Real> eigvals;
    stress.symmetricEigenvalues(eigvals);
    Real tmp = eigvals[2] - eigvals[0] + (eigvals[2] + eigvals[0] - 2*mean_stress)*sin_angle;
    Real dtmp = (eigvals[2] + eigvals[0] - 2*mean_stress)*dsin_angle;
    Real denom = std::sqrt(_small_smoother2 + 0.25*std::pow(tmp, 2));
    df_dintnl[0].assign(1, mean_stress*dsin_angle + 0.25*tmp*dtmp/denom - dcohesion(intnl[0])*cos_angle - cohesion(intnl[0])*dcos_angle);
  }
  else
  {
    // the edge-smoothed version
    Real aaa, bbb, ccc;
    abbo(sin3Lode, sin_angle, aaa, bbb, ccc);
    Real daaa, dbbb, dccc;
    dabbo(sin3Lode, sin_angle, daaa, dbbb, dccc);
    Real kk = aaa + bbb*sin3Lode + ccc*std::pow(sin3Lode, 2);
    Real dkk = (daaa + dbbb*sin3Lode + dccc*std::pow(sin3Lode, 2))*dsin_angle;
    Real sibar2 = stress.secondInvariant();
    Real denom = std::sqrt(_small_smoother2 + sibar2*std::pow(kk, 2));
    df_dintnl[0].assign(1, mean_stress*dsin_angle + sibar2*kk*dkk/denom - dcohesion(intnl[0])*cos_angle - cohesion(intnl[0])*dcos_angle);
  }
}

void
FiniteStrainCriticalStateMohrCoulomb::flowPotential(const RankTwoTensor & stress, const std::vector<Real> & intnl, std::vector<RankTwoTensor> & r)
{
  Real p = stress.trace()/3.0;
  Real sinpsi = std::sin(psi(intnl[0], stress));
  r.assign(1, df_dsig(stress, sinpsi));
}

void
FiniteStrainCriticalStateMohrCoulomb::dflowPotential_dstress(const RankTwoTensor & stress, const std::vector<Real> & intnl, std::vector<RankFourTensor> & dr_dstress)
{
  Real sin_angle = std::sin(psi(intnl[0], stress));
  Real mean_stress = stress.trace()/3.0;
  RankTwoTensor dmean_stress = stress.dtrace()/3.0;
  Real sin3Lode = stress.sin3Lode(_lode_cutoff, 0);
  if (std::abs(sin3Lode) <= _sin3tt)
  {
    // the non-edge-smoothed version
    std::vector<Real> eigvals;
    std::vector<RankTwoTensor> deigvals;
    std::vector<RankFourTensor> d2eigvals;
    stress.dsymmetricEigenvalues(eigvals, deigvals);
    stress.d2symmetricEigenvalues(d2eigvals);

    Real tmp = eigvals[2] - eigvals[0] + (eigvals[2] + eigvals[0] - 2*mean_stress)*sin_angle;
    RankTwoTensor dtmp = deigvals[2] - deigvals[0] + (deigvals[2] + deigvals[0] - 2*dmean_stress)*sin_angle;
    Real denom = std::sqrt(_small_smoother2 + 0.25*std::pow(tmp, 2));

    dr_dstress.assign(1, 0.25*tmp*(d2eigvals[2] - d2eigvals[0] + (d2eigvals[2] + d2eigvals[0])*sin_angle)/denom);
    Real pre = (0.25 - std::pow(0.25*tmp/denom, 2))/denom;
    for (unsigned i = 0; i < 3; ++i)
      for (unsigned j = 0; j < 3; ++j)
        for (unsigned k = 0; k < 3; ++k)
          for (unsigned l = 0; l < 3; ++l)
            dr_dstress[0](i, j, k, l) += pre*dtmp(i, j)*dtmp(k, l);
  }
  else
  {
    // the edge-smoothed version
    Real aaa, bbb, ccc;
    abbo(sin3Lode, sin_angle, aaa, bbb, ccc);
    RankTwoTensor dsin3Lode = stress.dsin3Lode(_lode_cutoff);
    Real kk = aaa + bbb*sin3Lode + ccc*std::pow(sin3Lode, 2);
    RankTwoTensor dkk = (bbb + 2*ccc*sin3Lode)*dsin3Lode;
    RankFourTensor d2kk = (bbb + 2*ccc*sin3Lode)*stress.d2sin3Lode(_lode_cutoff);
    for (unsigned i = 0; i < 3; ++i)
      for (unsigned j = 0; j < 3; ++j)
        for (unsigned k = 0; k < 3; ++k)
          for (unsigned l = 0; l < 3; ++l)
            d2kk(i, j, k, l) += 2*ccc*dsin3Lode(i, j)*dsin3Lode(k, l);

    Real sibar2 = stress.secondInvariant();
    RankTwoTensor dsibar2 = stress.dsecondInvariant();
    RankFourTensor d2sibar2 = stress.d2secondInvariant();

    Real denom = std::sqrt(_small_smoother2 + sibar2*std::pow(kk, 2));
    dr_dstress.assign(1, (0.5*d2sibar2*std::pow(kk, 2) + sibar2*kk*d2kk)/denom);
    for (unsigned i = 0; i < 3; ++i)
      for (unsigned j = 0; j < 3; ++j)
        for (unsigned k = 0; k < 3; ++k)
          for (unsigned l = 0; l < 3; ++l)
          {
            dr_dstress[0](i, j, k, l) += (dsibar2(i, j)*dkk(k, l)*kk + dkk(i, j)*dsibar2(k, l)*kk + sibar2*dkk(i, j)*dkk(k, l))/denom;
            dr_dstress[0](i, j, k, l) -= (0.5*dsibar2(i, j)*std::pow(kk, 2) + sibar2*kk*dkk(i, j))*(0.5*dsibar2(k, l)*std::pow(kk, 2) + sibar2*kk*dkk(k, l))/std::pow(denom, 3);
          }
  }

}

void
FiniteStrainCriticalStateMohrCoulomb::dflowPotential_dintnl(const RankTwoTensor & stress, const std::vector<Real> & intnl, std::vector<std::vector<RankTwoTensor> > & dr_dintnl)
{
  dr_dintnl.resize(1);

  Real sin_angle = std::sin(psi(intnl[0], stress));
  Real dsin_angle = std::cos(psi(intnl[0], stress))*dpsi(intnl[0]);

  Real mean_stress = stress.trace()/3.0;
  RankTwoTensor dmean_stress = stress.dtrace()/3.0;
  Real sin3Lode = stress.sin3Lode(_lode_cutoff, 0);

  if (std::abs(sin3Lode) <= _sin3tt)
  {
    // the non-edge-smoothed version
    std::vector<Real> eigvals;
    std::vector<RankTwoTensor> deigvals;
    stress.dsymmetricEigenvalues(eigvals, deigvals);
    Real tmp = eigvals[2] - eigvals[0] + (eigvals[2] + eigvals[0] - 2*mean_stress)*sin_angle;
    Real dtmp_dintnl = (eigvals[2] + eigvals[0] - 2*mean_stress)*dsin_angle;
    RankTwoTensor dtmp_dstress = deigvals[2] - deigvals[0] + (deigvals[2] + deigvals[0] - 2*dmean_stress)*sin_angle;
    RankTwoTensor d2tmp_dstress_dintnl = (deigvals[2] + deigvals[0] - 2*dmean_stress)*dsin_angle;
    Real denom = std::sqrt(_small_smoother2 + 0.25*std::pow(tmp, 2));
    dr_dintnl[0].assign(1, dmean_stress*dsin_angle + 0.25*dtmp_dintnl*dtmp_dstress/denom + 0.25*tmp*d2tmp_dstress_dintnl/denom - 0.25*tmp*dtmp_dstress*0.25*tmp*dtmp_dintnl/std::pow(denom, 3));
  }
  else
  {
    // the edge-smoothed version
    Real aaa, bbb, ccc;
    abbo(sin3Lode, sin_angle, aaa, bbb, ccc);
    Real kk = aaa + bbb*sin3Lode + ccc*std::pow(sin3Lode, 2);

    Real daaa, dbbb, dccc;
    dabbo(sin3Lode, sin_angle, daaa, dbbb, dccc);
    Real dkk_dintnl = (daaa + dbbb*sin3Lode + dccc*std::pow(sin3Lode, 2))*dsin_angle;
    RankTwoTensor dkk_dstress = (bbb + 2*ccc*sin3Lode)*stress.dsin3Lode(_lode_cutoff);
    RankTwoTensor d2kk_dstress_dintnl = (dbbb + 2*dccc*sin3Lode)*stress.dsin3Lode(_lode_cutoff)*dsin_angle;

    Real sibar2 = stress.secondInvariant();
    RankTwoTensor dsibar2 = stress.dsecondInvariant();
    Real denom = std::sqrt(_small_smoother2 + sibar2*std::pow(kk, 2));

    dr_dintnl[0].assign(1, dmean_stress*dsin_angle + (dsibar2*kk*dkk_dintnl + sibar2*dkk_dintnl*dkk_dstress + sibar2*kk*d2kk_dstress_dintnl)/denom - (0.5*dsibar2*std::pow(kk, 2) + sibar2*kk*dkk_dstress)*sibar2*kk*dkk_dintnl/std::pow(denom, 3));
  }
}

void
FiniteStrainCriticalStateMohrCoulomb::hardPotential(const RankTwoTensor & /*stress*/, const std::vector<Real> & /*intnl*/, std::vector<std::vector<Real> > & h)
{
  h.resize(1);
  h[0].assign(1, -1.0);
}

void
FiniteStrainCriticalStateMohrCoulomb::dhardPotential_dstress(const RankTwoTensor & /*stress*/, const std::vector<Real> & /*intnl*/, std::vector<std::vector<RankTwoTensor> > & dh_dstress)
{
  dh_dstress.resize(1);
  dh_dstress[0].assign(1, RankTwoTensor());
}

void
FiniteStrainCriticalStateMohrCoulomb::dhardPotential_dintnl(const RankTwoTensor & /*stress*/, const std::vector<Real> & /*intnl*/, std::vector<std::vector<std::vector<Real> > > & dh_dintnl)
{
  dh_dintnl.resize(1);
  dh_dintnl[0].resize(1);
  dh_dintnl[0][0].assign(1, 0.0);
}


Real
FiniteStrainCriticalStateMohrCoulomb::cohesion(const Real internal_param)
{
  return _cohesion_residual + (_cohesion - _cohesion_residual)*std::exp(-_cohesion_rate*internal_param);
}

Real
FiniteStrainCriticalStateMohrCoulomb::dcohesion(const Real internal_param)
{
  return -_cohesion_rate*(_cohesion - _cohesion_residual)*std::exp(-_cohesion_rate*internal_param);
}

Real
FiniteStrainCriticalStateMohrCoulomb::phi(const Real internal_param)
{
  return _phi_residual + (_phi - _phi_residual)*std::exp(-_phi_rate*internal_param);
}

Real
FiniteStrainCriticalStateMohrCoulomb::dphi(const Real internal_param)
{
  return -_phi_rate*(_phi - _phi_residual)*std::exp(-_phi_rate*internal_param);
}

Real
FiniteStrainCriticalStateMohrCoulomb::psi(const Real internal_param, const RankTwoTensor & stress)
{
  Real p = std::max(-1.0 * stress.trace()/3.0, 0.0);
  return _psi_0 * (_gamma - _lambda*std::pow(p/_p_ref, _xi) - _void_ratio[_qp]);
}

Real
FiniteStrainCriticalStateMohrCoulomb::dpsi(const Real internal_param)
{
  return 0;
}


void
FiniteStrainCriticalStateMohrCoulomb::abbo(const Real sin3lode, const Real sin_angle, Real & aaa, Real & bbb, Real & ccc)
{
  Real tmp1 = (sin3lode >= 0 ? _costt - sin_angle*_sintt/std::sqrt(3.0) : _costt + sin_angle*_sintt/std::sqrt(3.0));
  Real tmp2 = (sin3lode >= 0 ? _sintt + sin_angle*_costt/std::sqrt(3.0) : -_sintt + sin_angle*_costt/std::sqrt(3.0));

  ccc = -_cos3tt*tmp1;
  ccc += (sin3lode >= 0 ? -3*_sin3tt*tmp2 : 3*_sin3tt*tmp2);
  ccc /= 18*std::pow(_cos3tt, 3);

  bbb = (sin3lode >= 0 ? _sin6tt*tmp1 : -_sin6tt*tmp1);
  bbb -= 6*_cos6tt*tmp2;
  bbb /= 18*std::pow(_cos3tt, 3);

  aaa = (sin3lode >= 0 ? -sin_angle*_sintt/std::sqrt(3.0) - bbb*_sin3tt : sin_angle*_sintt/std::sqrt(3.0) + bbb*_sin3tt);
  aaa += -ccc*std::pow(_sin3tt, 2) + _costt;
}

void
FiniteStrainCriticalStateMohrCoulomb::dabbo(const Real sin3lode, const Real /*sin_angle*/, Real & daaa, Real & dbbb, Real & dccc)
{
  Real dtmp1 = (sin3lode >= 0 ? -_sintt/std::sqrt(3.0) : _sintt/std::sqrt(3.0));
  Real dtmp2 = _costt/std::sqrt(3.0);

  dccc = -_cos3tt*dtmp1;
  dccc += (sin3lode >= 0 ? -3*_sin3tt*dtmp2 : 3*_sin3tt*dtmp2);
  dccc /= 18*std::pow(_cos3tt, 3);

  dbbb = (sin3lode >= 0 ? _sin6tt*dtmp1 : -_sin6tt*dtmp1);
  dbbb -= 6*_cos6tt*dtmp2;
  dbbb /= 18*std::pow(_cos3tt, 3);

  daaa = (sin3lode >= 0 ? -_sintt/std::sqrt(3.0) - dbbb*_sin3tt : _sintt/std::sqrt(3.0) + dbbb*_sin3tt);
  daaa += -dccc*std::pow(_sin3tt, 2);
}


// set elasticity tensor from given stress and void ratio
RankFourTensor
FiniteStrainCriticalStateMohrCoulomb::get_elasticity_tensor(const RankTwoTensor & stress, const Real void_ratio)
{
  // set current mean effective stress
  Real p = -stress.trace()/3.0;
  p = std::max(p, 5.0);

  // set shear modulus and lambda
  Real shear_modulus = _c_g * (2.17 - void_ratio) * (2.17 - void_ratio) / (1.0 + void_ratio) * std::sqrt(_p_ref * p);
  Real lambda = 2.0 * shear_modulus * _poisson_ratio / (1.0 - 2 * _poisson_ratio);

  // set elasticity vector
  std::vector<Real> tmp;
  tmp.assign(2, 0);
  tmp[0] = lambda;
  tmp[1] = shear_modulus;

  RankFourTensor cijkl;
  cijkl.fillFromInputVector(tmp, (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"));

  return cijkl;
}
