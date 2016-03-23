#ifndef LONGTAILEDTITAPP_H
#define LONGTAILEDTITAPP_H

#include "MooseApp.h"

class LongtailedtitApp;

template<>
InputParameters validParams<LongtailedtitApp>();

class LongtailedtitApp : public MooseApp
{
public:
  LongtailedtitApp(InputParameters parameters);
  virtual ~LongtailedtitApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* LONGTAILEDTITAPP_H */
