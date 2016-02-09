

namespace jags {
  namespace covafillJAGS {

    class covafillModule : public Module {
    public:
      covafillModule();
      ~covafillModule();
    };

    covafillModule::covafillModule() : Module("covafillr") {
      insert(new covafillJAGS);
    }

    covafillModule::~covafillModule() {
      vector<Function*> const &fvec = functions();
      for (unsigned int i = 0; i < fvec.size(); ++i) {
	delete fvec[i];
      }
    }
  }
}

jags::covafillJAGS::covafillModule _covafillr_module;
