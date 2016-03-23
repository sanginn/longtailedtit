#include "LongtailedtitApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template<>
InputParameters validParams<LongtailedtitApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;

  return params;
}

LongtailedtitApp::LongtailedtitApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  LongtailedtitApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  LongtailedtitApp::associateSyntax(_syntax, _action_factory);
}

LongtailedtitApp::~LongtailedtitApp()
{
}

// External entry point for dynamic application loading
extern "C" void LongtailedtitApp__registerApps() { LongtailedtitApp::registerApps(); }
void
LongtailedtitApp::registerApps()
{
  registerApp(LongtailedtitApp);
}

// External entry point for dynamic object registration
extern "C" void LongtailedtitApp__registerObjects(Factory & factory) { LongtailedtitApp::registerObjects(factory); }
void
LongtailedtitApp::registerObjects(Factory & factory)
{
}

// External entry point for dynamic syntax association
extern "C" void LongtailedtitApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { LongtailedtitApp::associateSyntax(syntax, action_factory); }
void
LongtailedtitApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
