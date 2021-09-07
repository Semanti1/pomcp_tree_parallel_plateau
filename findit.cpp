#include "findit.h"
#include "utils.h"
#include <tuple>
using namespace std;
using namespace UTILS;

// NOTE: Very preliminary, currently, there isn't much "intelligence" (aka, we're not currently really 
// using the statistic that keeps track of the probability that the object is here--class variable for FINDIT_STATE)

FINDIT::FINDIT(int xsize, int ysize)
:   sensor_uncertainty(0.1),
    board_size(xsize)
{
    NumActions = 8; //(xsize * xsize) + 1;
    // 0: PICK_UP
    // 1 - xsize = coordinates
    NumObservations = 3;
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


STATE* FINDIT::Copy(const STATE& state) const {
    const FINDIT_STATE& fi_state = safe_cast<const FINDIT_STATE&>(state);
    FINDIT_STATE* new_state = MemoryPool.Allocate();
    *new_state = fi_state;
    
    return new_state;
}

void FINDIT::Validate(const STATE& state) const {
    const FINDIT_STATE& fi_state = safe_cast<const FINDIT_STATE&>(state);
    //assert(fi_state.object_position.X < board_size);
    //assert(fi_state.object_position.Y < board_size);
    //assert(fi_state)
}

STATE* FINDIT::CreateStartState() const {
    FINDIT_STATE* fi_state = MemoryPool.Allocate();
    //fi_state->positions.Resize(board_size, board_size);

    //fi_state->agent_position = COORD(0, 0);
    fi_state->agent_position = std::make_tuple(2, 0, 2);
    fi_state->orientation_r = 90;
    fi_state->orientation_o = Random(90);
    fi_state->foundit = false;
    fi_state->oriented = false;
    /*for (int i = 0; i < (board_size * board_size); i++) { 
        FINDIT_STATE::CELL& cell = fi_state->positions(i);
        cell.occupied = false;
        cell.visited = false;
        cell.likelihood_hit = 1.0;
        cell.likelihood_miss = 1.0;
        cell.probability_obj_here = (1.0 / (board_size * board_size));
    }*/

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

void FINDIT::FreeState(STATE* state) const {
    FINDIT_STATE* fi_state = safe_cast<FINDIT_STATE*>(state);
    MemoryPool.Free(fi_state);
}

/**
 * Convert coordinate to corresponding action
 */
/*int FINDIT::AgentPosToAction(COORD agent_pos) const {
    return 1 + agent_pos.X + (agent_pos.Y * board_size);
}*/

bool FINDIT::Step(STATE& state, int action, int& observation, double& reward) const {
    FINDIT_STATE& fi_state = safe_cast<FINDIT_STATE&>(state);
    reward = 0;
    observation = OBS_NONE;
    
    /*Sreshtaa's code */
    /*switch (action) {
        case PICK_UP: {
            if (fi_state.agent_position == fi_state.object_position) {
                reward = +100;
            } else {
                reward = -100;
            }

            return true;
	    }

        default: {
            // We'll probably want to change this to use the robot's input, but for now I encoded it in software using 
            // a 10% error rate (encoded by sensor_uncertainty)
            fi_state.agent_position = COORD((action - 1) % board_size, (action - 1) / board_size);

            FINDIT_STATE::CELL& visited_cell = fi_state.positions(action);
            if (visited_cell.visited) { // Discourages robot from continually touching the same cell
                reward = -3;
            } else {
                reward = -1;
                visited_cell.visited = true;
            }

            if (Bernoulli(sensor_uncertainty)) {
                observation = (fi_state.agent_position == fi_state.object_position) ? OBS_MISS : OBS_HIT;
            } else {
                observation = (fi_state.agent_position == fi_state.object_position) ? OBS_HIT : OBS_MISS;
            }

            // Do some stats thing here to take into account the uncertainty in the observation (such that having multiple
            // observations that report the same thing will increase the amount of confidence we have) 
            // Don't think this following if block is necessary (because I don't think recieved_hit_obs is necessary)
            if (observation == OBS_HIT) {
                visited_cell.likelihood_hit *= (1.0 - sensor_uncertainty);
                visited_cell.likelihood_miss *= sensor_uncertainty;
            } else {
                visited_cell.likelihood_hit *= sensor_uncertainty;
                visited_cell.likelihood_miss *= (1.0 - sensor_uncertainty);
            }
            
            // taken from rocksample
            double denom = (0.5 * visited_cell.likelihood_miss) + (0.5 * visited_cell.likelihood_hit);
		    visited_cell.probability_obj_here = (0.5 * visited_cell.likelihood_hit) / denom;

            break;
        }
    }*/

    /* Semanti's code*/
    std::tuple<int, int, int> r_pos = fi_state.agent_position;

    int rpos_x = get<0>(r_pos);
    int rpos_y = get<1>(r_pos);
    int rpos_z = get<2>(r_pos);
    switch (action) {
        

            

    case MOVE_EAST: {
        if (rpos_x + 1 < board_size)
        {

            fi_state.agent_position = std::make_tuple(rpos_x + 1, rpos_y, 1);
           // fi_state.agent_position = std::make_tuple(rpos_x, rpos_y, 1);
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
            //reward = reward - 1;
            break;
           // reward = reward - 1;
           // break;
        }
       // cout<< "Agent pos "<< fi_state.agent_position
        break; }

    case MOVE_WEST: {
        if (rpos_x - 1 >= 0)
        {
            fi_state.agent_position = std::make_tuple(rpos_x - 1, rpos_y, 1);
            //fi_state.agent_position = std::make_tuple(rpos_x, rpos_y, 1);
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
            
            break;
            //reward = reward - 1;
           // break;
        }
        break; }

    case MOVE_NORTH: {
        if (rpos_y + 1 < board_size)
        {
            fi_state.agent_position = std::make_tuple(rpos_x, rpos_y + 1, 1);
            //fi_state.agent_position = std::make_tuple(rpos_x, rpos_y, 1);
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
            //reward = reward - 1;
            break;
           // reward = reward - 1;
            //break;
        }
        break; }

        case MOVE_SOUTH: {
            if (rpos_y - 1 >= 0)
            {
                fi_state.agent_position = std::make_tuple(rpos_x, rpos_y - 1, 1);
                //fi_state.agent_position = std::make_tuple(rpos_x, rpos_y, 1);
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
                break;
                //reward = reward - 1;
                //break;
            }
            break; }
        /*case LOOK:
        {
            fi_state.agent_position = std::make_tuple(rpos_x, rpos_y, 1);
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
        reward = reward - 1;
        break; 
        }*/
        /*case FOUND: {
            //fi_state.foundit = true;
            if (fi_state.agent_position == fi_state.object_position)
            {
                fi_state.foundit = true;
                reward += 50;
            }
            else
                reward -= 50;
            //observation = GetObservation();
            break; }*/
        /*case ORIENT: {
            fi_state.orientation_r = Random(90);
            reward = reward - 1;
            break; }
        case ORIENT_DONE: {
            fi_state.oriented = true;
            if (abs(fi_state.orientation_o - fi_state.orientation_r) <= 5)
            {
                reward = reward + 50;
                return true;
            }
            else
                reward = reward - 50;
            break; }
        case DONE: {
            if ((fi_state.agent_position == fi_state.object_position) && abs(fi_state.orientation_o - fi_state.orientation_r) <= 5)
            {
                reward = reward + 100;
                return true;
            }
            else
                reward = reward - 100;
            break; }*/




    }

    return false;
}

// ############################################################################

/**
 * Honestly not so sure of what this is supposed to do
 */
bool FINDIT::LocalMove(STATE& state, const HISTORY& history, int stepObs, const STATUS& status) const {
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
void FINDIT::GenerateLegal(const STATE& state, const HISTORY& history, 
                            std::vector<int>& legal, const STATUS& status) const {
    
    const FINDIT_STATE& fi_state = safe_cast<const FINDIT_STATE&>(state);

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
        std::tuple<int,int,int> r_pos = fi_state.agent_position;
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
            legal.push_back(MOVE_SOUTH);

        //legal.push_back(FOUND);
       // legal.push_back(LOOK);
  //  }

       

}

// ############################################################################

void FINDIT::DisplayBeliefs(const BELIEF_STATE& beliefState, std::ostream& ostr) const {
    // Not implemented
    return;
}

void FINDIT::DisplayState(const STATE& state, std::ostream& ostr) const {
    // Not implemented
    const FINDIT_STATE& fi_state = safe_cast<const FINDIT_STATE&>(state);
    ostr << "Agent pos " <<get<0>(fi_state.agent_position) << ","<< get<1>(fi_state.agent_position)<< ","<< get<2>(fi_state.agent_position) << endl;
    ostr << "Object pos " << get<0>(fi_state.object_position) <<","<< get<1>(fi_state.object_position) << "," << get<2>(fi_state.object_position) <<endl;
    ostr << "Found it? " << fi_state.foundit << endl;
    ostr << "Agent orientation " << fi_state.orientation_r<< endl;
    ostr << "obj orientation " << fi_state.orientation_o << endl;
    ostr << "Oriented it? " << fi_state.oriented << endl;
   // cout << 

    //return;
}

void FINDIT::DisplayObservation(const STATE& state, int observation, std::ostream& ostr) const {
    // Not implemented
    switch (observation)
    {
    case OBS_NONE:
        break;
    case OBS_HIT:
        ostr << "Observed pos" << endl;
        break;
    case OBS_MISS:
        ostr << "Observed negative" << endl;
        break;
    }
  //  return;
}

void FINDIT::DisplayAction(int action, std::ostream& ostr) const {
    // Not implemented
    switch (action) 
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
        break; }
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
}
