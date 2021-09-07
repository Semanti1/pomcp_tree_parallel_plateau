#ifndef FINDIT_IMPROVED_H
#define FINDIT_IMPROVED_H
#include <tuple>
#include "simulator.h"
#include "grid.h"
#include "coord.h"
#include <list>

// NOTE: Very preliminary, there are definitely places for improvement! 
class FINDIT_IMPROVED_STATE : public STATE {
public: 
    struct CELL {
       // bool occupied;
        bool Visited;
        //double likelihood_hit; // used for calculation of probability_obj_here
        //double likelihood_miss; // used for calculation of probability_obj_here
        //double probability_obj_here;
    };

    GRID<CELL> positions;
    /*COORD object_position;
    COORD agent_position;*/
    std::tuple<int, int, int> object_position;
    std::tuple<int, int, int> agent_position;
    bool foundit;
    //int orientation_r;
    //int orientation_o;
    //bool oriented;

};


class FINDIT_IMPROVED : public SIMULATOR {
public:
    FINDIT_IMPROVED(int xsize, int ysize);
    virtual STATE* Copy(const STATE& state) const;
    virtual void Validate(const STATE& state) const;
    virtual STATE* CreateStartState() const;
    virtual void FreeState(STATE* state) const;
    virtual bool Step(STATE& state, int action, 
        int& observation, double& reward) const;

    virtual bool LocalMove(STATE& state, const HISTORY& history,
        int stepObs, const STATUS& status) const;
    void GenerateLegal(const STATE& state, const HISTORY& history,
        std::vector<int>& legal, const STATUS& status) const;
    // void GeneratePreferred(const STATE& state, const HISTORY& history,
    //     std::vector<int>& preferred, const STATUS& status) const;

    virtual void DisplayBeliefs(const BELIEF_STATE& beliefState, 
        std::ostream& ostr) const;
    virtual void DisplayState(const STATE& state, std::ostream& ostr) const;
    virtual void DisplayObservation(const STATE& state, int observation, std::ostream& ostr) const;
    virtual void DisplayAction(int action, std::ostream& ostr) const;

protected: 
    int board_size;
    double sensor_uncertainty;
    
    // Observations
    enum {
        OBS_HIT = 0, 
        OBS_MISS = 1, 
     //   OBS_NONE = 2
    };

    // Extra actions (apart from the coordinate ones)
   /* enum {
        MOVE_EAST = 0,
        MOVE_WEST = 1,
        MOVE_NORTH = 2,
        MOVE_SOUTH = 3,*/
      //  LOOK = 4,
        //FOUND = 4,
        //ORIENT = 5,
        //ORIENT_DONE = 6,
        //DONE = 7
        //PICK_UP = 0
  //  };
    
private: 
    mutable MEMORY_POOL<FINDIT_IMPROVED_STATE> MemoryPool;
    //int AgentPosToAction(COORD agent_pos) const;

    // COORD object_position;
};

#endif
