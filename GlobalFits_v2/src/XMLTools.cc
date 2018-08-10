#include "MiniBooNE.h"
#include "LSND_loglikelihood.h"
#include "Gallium.h"
#include "NOMAD.h"
#include "CDHS.h"
#include "CCFR.h"
#include "MiniBooNE_dis.h"
#include "NuMI.h"
#include "KARMEN.h"
#include "XMLTools.h"

int FitReader::Load(std::string xml){

  // Load XML file
  XMLDocument doc;
  if (doc.LoadFile(xml.c_str())){
    std::cout << "Couldn't load XML file! I quit." << std::endl;
    return 1;
  }
  XMLHandle hDoc(&doc);

  // We'll have an oscillator and several datasets
  XMLElement *pData, *pOsc;
  pOsc = doc.FirstChildElement("oscillator");
  pData = doc.FirstChildElement("dataset");

  std::string dset;

  if (!pData){
    std::cout << "No datasets in config. Outta here." << std::endl;
    return 1;
  }
  else while(pData){
    dset = pData->Attribute("name");
    if(stoi(pData->Attribute("use"))){
      if(dset == "MBnu"){
        myDataSets.push_back(new MiniBooNE(false));
        std::cout << "Using MiniBooNE Nu dataset" << std::endl;
      }
      else if(dset == "MBnubar"){
        myDataSets.push_back(new MiniBooNE(true));
        std::cout << "Using MiniBooNE Nubar dataset" << std::endl;
      }
      else if(dset == "NuMI"){
        myDataSets.push_back(new NuMI);
        std::cout << "Using MiniBooNE NuMI Dataset" << std::endl;
      }
      else if(dset == "LSND_loglikelihood"){
        myDataSets.push_back(new LSND_loglikelihood);
        std::cout << "Using LSND (log likelihood mode)" << std::endl;
      }
      else if(dset == "CCFR"){
        myDataSets.push_back(new CCFR);
        std::cout << "Using CCFR" << std::endl;
      }
      else if(dset == "Gallium"){
        myDataSets.push_back(new Gallium);
        std::cout << "Using Gallium" << std::endl;
      }
      else if(dset == "MBnu_dis"){
        myDataSets.push_back(new MiniBooNE_dis(false));
        std::cout << "Using MiniBooNE Nu Disappearance dataset" << std::endl;
      }
      else if(dset == "MBnubar_dis"){
        myDataSets.push_back(new MiniBooNE_dis(true));
        std::cout << "Using MiniBooNE Nubar Disappearance dataset" << std::endl;
      }
      else if(dset == "NOMAD"){
        myDataSets.push_back(new NOMAD);
        std::cout << "Using NOMAD dataset" << std::endl;
      }
      else if(dset == "CDHS"){
        myDataSets.push_back(new CDHS);
        std::cout << "Using CDHS dataset" << std::endl;
      }
      else if(dset == "KARMEN"){
        myDataSets.push_back(new KARMEN);
        std::cout << "Using KARMEN dataset" << std::endl;
      }
      else
        std::cout << "Dataset not implemented yet!" << std::endl;
    }
    pData = pData->NextSiblingElement("dataset");
  }
  if(myDataSets.size()==0){
    std::cout << "No valid datasets requested." << std::endl;
    return 1;
  }

  if (!pOsc) {
    std::cout << "No oscillator configuration. I'm outta here."  << std::endl;
    return 1;
  }
  else{
    int nsteriles, cpcons, rndseed, grdpts;
    float umax, usqmax, stepsize, temperature;
    nsteriles = stoi(pOsc->Attribute("nsteriles"));
    cpcons = stoi(pOsc->Attribute("cpcons"));
    rndseed = stoi(pOsc->Attribute("rndseed"));
    grdpts = stoi(pOsc->Attribute("gridpts"));
    umax = stof(pOsc->Attribute("umax"));
    usqmax = stof(pOsc->Attribute("usqmax"));
    stepsize = stof(pOsc->Attribute("stepsize"));
    temperature = stof(pOsc->Attribute("temperature"));

    myOscillator = Oscillator(.01f,100.f,0.f,umax,usqmax,stepsize,temperature,nsteriles,grdpts,cpcons,rndseed);
  }

  return 0;
}

int ProcessReader::Load(std::string xml){

  // Load XML file
  XMLDocument doc;
  if (doc.LoadFile(xml.c_str())){
    std::cout << "Couldn't load XML file! I quit." << std::endl;
    return 1;
  }
  XMLHandle hDoc(&doc);

  // We'll have an oscillator and several datasets
  XMLElement *pData, *pProc;
  pData = doc.FirstChildElement("dataset");
  pProc = doc.FirstChildElement("procopts");

  std::string dset,dloc;

  if (!pData){
    std::cout << "No datasets in config. Outta here." << std::endl;
    return 1;
  }
  else while(pData){
    dset = pData->Attribute("name");
    if(stoi(pData->Attribute("use"))){
      raster = stoi(pData->Attribute("raster"));
      dloc = pData->Attribute("loc");
      std::cout << "Adding tree " << dset << " from " << dloc << "." << std::endl;
      data_names.push_back(dset);
      data_files.push_back(new TFile(dloc.c_str(),"READ"));
      data_trees.push_back((TTree*)data_files.back()->Get(dset.c_str()));
    }

    pData = pData->NextSiblingElement("dataset");
  }
  if(data_names.size()==0){
    std::cout << "No valid datasets requested." << std::endl;
    return 1;
  }

  if (!pProc){
    std::cout << "No process options. Fuck this." << std::endl;
    return 1;
  }
  else{
    tag = pProc->Attribute("tag");
    gridpts_dm2 = atoi(pProc->Attribute("gridpts_dm2"));
    gridpts_sin22th = atoi(pProc->Attribute("gridpts_sin22th"));
  }
  return 0;
}
