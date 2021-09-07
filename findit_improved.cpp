#include "findit_improved.h"
#include "utils.h"
#include <tuple>
using namespace std;
using namespace UTILS;

// NOTE: Very preliminary, currently, there isn't much "intelligence" (aka, we're not currently really 
// using the statistic that keeps track of the probability that the object is here--class variable for FINDIT_STATE)

FINDIT_IMPROVED::FINDIT_IMPROVED(int xsize, int ysize)
:   sensor_uncertainty(0.1),
    board_size(xsize)
{
    NumActions =(xsize * xsize);
    // 0: PICK_UP
    // 1 - xsize = coordinates
    NumObservations = 2;
    // 0: OBS_HIT
    // 1: OBS_MISS
    // 2: OBS_NONE
    RewardRange = 100;
    Discount = 0.95;

    // RandomSeed((unsigned) time(0)); 
    // int obj_pos_x = Random(board_size);
    // int obj_pos_y = Random(board_size);
    // object_position = COORD(obj_pos_x, obj_pos_y);
}


STATE* FINDIT_IMPROVED::Copy(const STATE& state) const {
    const FINDIT_IMPROVED_STATE& fi_state = safe_cast<const FINDIT_IMPROVED_STATE&>(state);
    FINDIT_IMPROVED_STATE* new_state = MemoryPool.Allocate();
    *new_state = fi_state;
    
    return new_state;
}

void FINDIT_IMPROVED::Validate(const STATE& state) const {
    const FINDIT_IMPROVED_STATE& fi_state = safe_cast<const FINDIT_IMPROVED_STATE&>(state);
    //assert(fi_state.object_position.X < board_size);
    //assert(fi_state.object_position.Y < board_size);
    //assert(fi_state)
}

STATE* FINDIT_IMPROVED::CreateStartState() const {
    FINDIT_IMPROVED_STATE* fi_state = MemoryPool.Allocate();
    fi_state->positions.Resize(board_size, board_size);

    //fi_state->agent_position = COORD(0, 0);
    fi_state->agent_position = std::make_tuple(2, 0, 2);
    //fi_state->orientation_r = 90;
    //fi_state->orientation_o = Random(90);
    fi_state->foundit = false;
    //fi_state->oriented = false;
    for (int i = 0; i < (board_size * board_size); i++) { 
        FINDIT_IMPROVED_STATE::CELL& cell = fi_state->positions(i);
      //  cell.occupied = false;
        cell.Visited = false;
        //cell.likelihood_hit = 1.0;
        //cell.likelihood_miss = 1.0;
        //cell.probability_obj_here = (1.0 / (board_size * board_size));
    }

    // I honestly am wondering if object position should be a private class variable instead of a state variable. 
    // In Rocksample, the rock positions are a private class variable, but in battleship, the ships are a state variable
    // It might be worth trying to move the next four lines to FINDIT::FINDIT(int xsize, int ysize)
    int obj_pos_x = Random(board_size);
    int obj_pos_y = Random(board_size);
    fi_state->object_position = std::make_tuple(obj_pos_x, obj_pos_y, 1);
    //fi_state->object_position = COORD(obj_pos_x, obj_pos_y); // if moving, change just to object_position and uncomment priv variable in header file
    // cout << "X = " << fi_state->object_position.X << " Y = " << fi_state->object_position.Y << endl;
   // fi_state->positions(fi_state->object_position).occupied = true;

    return fi_state;
}

void FINDIT_IMPROVED::FreeState(STATE* state) const {
    FINDIT_IMPROVED_STATE* fi_state = safe_cast<FINDIT_IMPROVED_STATE*>(state);
    MemoryPool.Free(fi_state);
}

/**
 * Convert coordinate to corresponding action
 */
/*int FINDIT::AgentPosToAction(COORD agent_pos) const {
    return 1 + agent_pos.X + (agent_pos.Y * board_size);
}*/

bool FINDIT_IMPROVED::Step(STATE& state, int action, int& observation, double& reward) const {
    FINDIT_IMPROVED_STATE& fi_state = safe_cast<FINDIT_IMPROVED_STATE&>(state);
    reward = 0;
    observation = 1;
    
    
    std::tuple<int, int, int> r_pos = fi_state.agent_position;

    int rpos_x = get<0>(r_pos);
    int rpos_y = get<1>(r_pos);
    int rpos_z = get<2>(r_pos);

    COORD actionPos = fi_state.positions.Coord(action);
   // cout << "Action POs" << actionPos;
    FINDIT_IMPROVED_STATE::CELL& cell = fi_state.positions(actionPos);
    if (cell.Visited)
    {
        reward = -10;
        observation = 1;
        //DisplayState(state, cout);
        //DisplayAction(action, cout);
      //  assert(false);
    }
    else
    {
        fi_state.agent_position = std::make_tuple(actionPos.X, actionPos.Y, 1);
        cell.Visited = true;
        double thresh = 0;
        if (fi_state.agent_position == fi_state.object_position)
            thresh = 1 - sensor_uncertainty;
        else
            thresh = sensor_uncertainty;

        double val = static_cast <double> (rand()) / static_cast <float> (RAND_MAX);
        if (val < thresh)
            observation = OBS_HIT;
        else
            observation = OBS_MISS;
        if (observation == 0)
        {
            if (fi_state.agent_position == fi_state.object_position)
            {
                fi_state.foundit = true;
                reward = reward + 100;
                return true;
            }
            else
            {
                reward = reward - 100;
            }
        }
        else
        {



            reward = reward - 1;
        }
    }
  

    return false;
}

// ############################################################################

/**
 * Honestly not so sure of what this is supposed to do
 */
bool FINDIT_IMPROVED::LocalMove(STATE& state, const HISTORY& history, int stepObs, const STATUS& status) const {
    return false;
    
    // Try and change object posiiton for particle reinvigoration
    // FINDIT_STATE& fi_state = safe_cast<FINDIT_STATE&>(state);
    // COORD current_agent_pos = fi_state.agent_position;
    // double probability_of_obj = fi_state.positions(current_agent_pos).probability_obj_here;

    // // Don't want to change object position if we're pretty sure object could be here 
    // if (probability_of_obj >= 0.75) { 
    //     return false;
    // }

    // COORD new_obj_pos = COORD(Random(board_size), Random(board_size));

    // if (fi_state.positions(new_obj_pos).visited) {
    //     return false;
    // } else {
    //     fi_state.object_position = new_obj_pos;
    //     return true;
    // }
}

/**
 * Strictly, what actions are we *allowed* to perform
 */
void FINDIT_IMPROVED::GenerateLegal(const STATE& state, const HISTORY& history,
                            std::vector<int>& legal, const STATUS& status) const {
    
    const FINDIT_IMPROVED_STATE& fi_state = safe_cast<const FINDIT_IMPROVED_STATE&>(state);
    for (int a = 0; a < NumActions; ++a)
        if (!fi_state.positions(a).Visited)
            legal.push_back(a);
    // COORD current_agent_pos = fi_state.agent_position;
    // double probability_of_obj = fi_state.positions(current_agent_pos).probability_obj_here;

    // We can always choose to pick up
    /*legal.push_back(PICK_UP);

    // We can pick any square in our grid to touch 
    for (int i = 1; i < (board_size + 1); i++) {
        legal.push_back(i);
    }*/
    

  /*  if (fi_state.oriented)
        legal.push_back(DONE);
    else if (fi_state.foundit)
    {
        legal.push_back(ORIENT);
        legal.push_back(ORIENT_DONE);
    }
    else {*/
      /*  std::tuple<int,int,int> r_pos = fi_state.agent_position;
        int rpos_x = std::get<0>(r_pos);
        int rpos_y = std::get<1>(r_pos);
        int rpos_z = std::get<2>(r_pos);

        if (rpos_x + 1 < board_size)
            legal.push_back(MOVE_EAST);

        if (rpos_y + 1 < board_size)
            legal.push_back(MOVE_NORTH);

        if (rpos_x - 1 >= 0)
            legal.push_back(MOVE_WEST);

        if (rpos_y - 1 >= 0)
            legal.push_back(MOVE_SOUTH);*/

        //legal.push_back(FOUND);
       // legal.push_back(LOOK);
  //  }

       

}

// ############################################################################

void FINDIT_IMPROVED::DisplayBeliefs(const BELIEF_STATE& beliefState, std::ostream& ostr) const {
    // Not implemented
    return;
}

void FINDIT_IMPROVED::DisplayState(const STATE& state, std::ostream& ostr) const {
    // Not implemented
    const FINDIT_IMPROVED_STATE& fi_state = safe_cast<const FINDIT_IMPROVED_STATE&>(state);
    ostr << "Agent pos " <<get<0>(fi_state.agent_position) << ","<< get<1>(fi_state.agent_position)<< ","<< get<2>(fi_state.agent_position) << endl;
    ostr << "Object pos " << get<0>(fi_state.object_position) <<","<< get<1>(fi_state.object_position) << "," << get<2>(fi_state.object_position) <<endl;
    ostr << "Found it? " << fi_state.foundit << endl;
   // ostr << "Agent orientation " << fi_state.orientation_r<< endl;
   // ostr << "obj orientation " << fi_state.orientation_o << endl;
   // ostr << "Oriented it? " << fi_state.oriented << endl;
   // cout << 

    //return;
}

void FINDIT_IMPROVED::DisplayObservation(const STATE& state, int observation, std::ostream& ostr) const {
    // Not implemented
    switch (observation)
    {
    
    case OBS_HIT:
        ostr << "Observed pos" << endl;
        break;
    case OBS_MISS:
        ostr << "Observed negative" << endl;
        break;
    }
  //  return;
}

void FINDIT_IMPROVED::DisplayAction(int action, std::ostream& ostr) const {
    // Not implemented

    COORD actionPos = COORD(action % board_size, action / board_size);
    ostr << "Action "<<actionPos << endl;
   /* ostr << endl << "  ";
    for (int x = 0; x < XSize; x++)
        ostr << setw(1) << x << ' ';
    ostr << "  " << endl;
    for (int y = YSize - 1; y >= 0; y--)
    {
        ostr << setw(1) << y << ' ';
        for (int x = 0; x < XSize; x++)
        {
            char c = ' ';
            if (actionPos == COORD(x, y))
                c = '@';
            ostr << c << ' ';
        }
        ostr << setw(1) << y << endl;
    }
    ostr << "  ";
    for (int x = 0; x < XSize; x++)
        ostr << setw(1) << x << ' ';
    ostr << "  " << endl;*/
   /* switch (action) 
    {
    case MOVE_EAST :
    {ostr << "Move East" << endl;
    break; }
    case MOVE_WEST: {
        ostr << "Move West" << endl;
        break; }
    case MOVE_NORTH: {
        ostr << "Move North" << endl;
        break; }
    case MOVE_SOUTH: {
        ostr << "Move South" << endl;
        break; }*/
    /*case LOOK: {
        ostr << "Look" << endl;
        break; }
    case FOUND: {
        ostr << "Found" << endl;
        break; }
    case  ORIENT: {
        ostr << "Orient" << endl;
        break; }
    case ORIENT_DONE: {
        ostr << "Orient Done" << endl;
        break; }
    case DONE: {
        ostr << "Done" << endl;
        break; }*/
    }
   // return;
//}
