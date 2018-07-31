#include "tinyxml2.h"
#include "datasets.h"
#include "OscTools.h"

using namespace tinyxml2;

class FitReader{
  public:
    FitReader(){};

    int Load(std::string xml);
    dataset* GetDataset(int ind) { return myDataSets[ind]; };
    int GetNDatasets()  { return myDataSets.size(); };
    Oscillator GetOscillator() { return  myOscillator; };

  private:
    std::vector < dataset* > myDataSets;
    Oscillator myOscillator;
};

class ProcessReader{
  public:
    ProcessReader(){};

    int Load(std::string xml);
    TTree* GetTree(int ind) { return data_trees[ind]; };
    int GetNTrees() { return data_trees.size(); };

    int gridpts_dm2, gridpts_sin22th;
    std::string tag;

  private:
    std::vector < std::string > data_names;
    std::vector < TTree* > data_trees;
    std::vector < TFile* > data_files;

};
