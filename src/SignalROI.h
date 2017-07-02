#ifndef WIRECELLSIGPROC_SignalROI_h
#define WIRECELLSIGPROC_SignalROI_h

#include "WireCellUtil/Waveform.h"

#include <iostream>
#include <vector>
#include <list>
#include <map>


namespace WireCell{
  namespace SigProc{
    class SignalROI{
    public:
      SignalROI(int plane, int chid, int start_bin, int end_bin, const Waveform::realseq_t& signal);
      SignalROI(SignalROI *roi);
      ~SignalROI();
      int get_start_bin(){return start_bin;}
      int get_end_bin(){return end_bin;}
      int get_chid(){return chid;}
      int get_plane(){return plane;}
      std::vector<float>& get_contents(){return contents;}
      std::vector<std::pair<int,int>> get_above_threshold(float th);
      
      bool overlap(SignalROI *roi);
      bool overlap(SignalROI *roi1, float th, float th1);
      
    private:
      int plane;
      int chid;
      int start_bin;
      int end_bin;
 
      
      std::vector<float> contents;
    };
    
    typedef std::list<SignalROI*>SignalROIList;
    typedef std::vector<SignalROI*> SignalROISelection; 
    typedef std::vector<SignalROISelection> SignalROIChSelection;
    typedef std::vector<SignalROIList> SignalROIChList;
    typedef std::map<SignalROI*, SignalROISelection> SignalROIMap;
  }
}

#endif
