#include "beliefstate.h"
#include "simulator.h"
#include "utils.h"

using namespace UTILS;

BELIEF_STATE::BELIEF_STATE()
{
   // std::mutex m;
   // std::lock_guard<std::mutex> lock(m);
    Samples.clear();
}

void BELIEF_STATE::Free(const SIMULATOR& simulator)
{
    //std::lock_guard<std::mutex> lock(m);
    for (std::vector<STATE*>::iterator i_state = Samples.begin();
            i_state != Samples.end(); ++i_state)
    {
        simulator.FreeState(*i_state);
    }
    Samples.clear();
}

STATE* BELIEF_STATE::CreateSample(const SIMULATOR& simulator, std::mutex* m) const
{
    //std::mutex m;
    std::lock_guard<std::mutex> lock(*m);
    int index = Random(Samples.size());
    STATE* retSt = simulator.Copy(*Samples[index]);
    //return simulator.Copy(*Samples[index]);
    return retSt;
}

void BELIEF_STATE::AddSamplem(STATE* state, std::mutex* m)
{
    //std::mutex m;
    std::lock_guard<std::mutex> lock(*m);
    Samples.push_back(state);
}
void BELIEF_STATE::AddSample(STATE* state)
{
    //std::mutex m;
    //std::lock_guard<std::mutex> lock(*m);
    Samples.push_back(state);
}

void BELIEF_STATE::Copy(const BELIEF_STATE& beliefs, const SIMULATOR& simulator)
{
    //std::mutex m;
    //std::lock_guard<std::mutex> lock(m);
    for (std::vector<STATE*>::const_iterator i_state = beliefs.Samples.begin();
        i_state != beliefs.Samples.end(); ++i_state)
    {
        //STATE* nst = 
        AddSample(simulator.Copy(**i_state));
    }
}

void BELIEF_STATE::Move(BELIEF_STATE& beliefs)
{
    //std::mutex m;
   // std::lock_guard<std::mutex> lock(mmove);
    for (std::vector<STATE*>::const_iterator i_state = beliefs.Samples.begin();
        i_state != beliefs.Samples.end(); ++i_state)
    {
        AddSample(*i_state);
    }
    beliefs.Samples.clear();
}
