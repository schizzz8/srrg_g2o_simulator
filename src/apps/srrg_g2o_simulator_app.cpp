#include <iostream>
#include <srrg_g2o_simulator/simulator.h>
#include <srrg_g2o_simulator_viewer/simulator_viewer.h>

#include <qapplication.h>
#include <qglviewer.h>

using namespace srrg_core;
using namespace srrg_matchable;
using namespace srrg_g2o_simulator;

int main(int argc, char **argv){

  float resolution = 1.0f;
  int d=6;
  int num_poses=10;

  int c=1;
  while (c<argc) {
    if(!strcmp(argv[c],"-d")){
      c++;
      d=std::atoi(argv[c]);
    } else if(!strcmp(argv[c],"-n")){
      c++;
      num_poses=std::atoi(argv[c]);
    } else if(!strcmp(argv[c],"-r")){
      c++;
      resolution=std::atof(argv[c]);
    }
    c++;
  }

  Simulator s;
  s.setResolution(resolution);
  s.setDimension(d);
  s.setNumPoses(num_poses);

  s.populateSet();

  s.sbragation();
  MatchablePtrSet set = s.matchables();

  Scene* scene = new Scene();
  for(MatchablePtrSet::iterator it = set.begin(); it != set.end(); ++it){
    MatchablePtr m = *it;
    scene->push_back(m);
  }

  std::cerr << "Scene has " << scene->size() << " matchables after sbraco" << std::endl;

  QApplication app(argc, argv);
  SimulatorViewer viewer;
  viewer.setScene(scene);
  viewer.show();
  app.exec();

  return 0;
}
