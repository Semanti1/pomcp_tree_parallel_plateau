#ifndef HISTORY_H
#define HISTORY_H

#include <vector>
#include <ostream>
#include <assert.h>
#include<mutex>
using namespace std;
#include <iostream>
//#include "mcts.h"
//extern std::mutex hmtx;

class HISTORY
{
public:

    struct ENTRY
    {
        ENTRY() { }

        ENTRY(int action, int obs)
        :   Action(action), Observation(obs)
        { }
        
        int Action;
        int Observation;
    };
    
    bool operator==(const HISTORY& history) const
    {
        //std::lock_guard<std::mutex> lock(m);
        //std::lock_guard<std::mutex> lock(hmtx);
        if (history.History.size() != History.size())
            return false;
        for (int i = 0; i < History.size(); ++i)
            if (history.History[i].Action != History[i].Action
             || history.History[i].Observation != History[i].Observation)
                return false;
        return true;
    }
    
    void Add(int action, int obs = -1) 
    { 
        //std::lock_guard<std::mutex> lock(hmtx);
        History.push_back(ENTRY(action, obs));
    }
    
    void Pop()
    {
       // std::lock_guard<std::mutex> lock(m);
       // std::lock_guard<std::mutex> lock(hmtx);
        History.pop_back();
    }
    
    void Truncate(int t)
    {
      //  std::lock_guard<std::mutex> lock(m);
       // std::lock_guard<std::mutex> lock(hmtx);
        History.resize(t);
    }
    
    void Clear() 
    { 
        //std::lock_guard<std::mutex> lock(m);
       // std::lock_guard<std::mutex> lock(hmtx);
        History.clear(); 
    }
    
    int Size() const
    {
       // std::lock_guard<std::mutex> lock(m);
      //  std::lock_guard<std::mutex> lock(hmtx);
        return History.size();
    }
    
    ENTRY& operator[](int t)
    {
       // std::lock_guard<std::mutex> lock(m);
      //  std::lock_guard<std::mutex> lock(hmtx);
       // cout << "T value range from "<< t << " to " << History.size() << endl;
        //#pragma omp atomic read
        assert(t >= 0 && t < History.size());
        return History[t];
    }

    const ENTRY& operator[](int t) const
    {
       // std::lock_guard<std::mutex> lock(m);
       // std::lock_guard<std::mutex> lock(hmtx);
        assert(t >= 0 && t < History.size());
        return History[t];
    }

    ENTRY& Back()
    {
      //  std::lock_guard<std::mutex> lock(m);
      //  std::lock_guard<std::mutex> lock(hmtx);
        assert(History.size() > 0);
        return History.back();
    }

    const ENTRY& Back() const
    {
       // std::lock_guard<std::mutex> lock(m);
      //  std::lock_guard<std::mutex> lock(hmtx);
        assert(History.size() > 0);
        return History.back();
    }

    void Display(std::ostream& ostr) const
    {
        //std::lock_guard<std::mutex> lock(m);
      //  std::lock_guard<std::mutex> lock(hmtx);
        for (int t = 0; t < History.size(); ++t)
        {
            ostr << "a=" << History[t].Action <<  " ";
            if (History[t].Observation >= 0)
                ostr << "o=" << History[t].Observation << " ";
        }
    }


private:

    std::vector<ENTRY> History;
    // std::mutex m;
};

#endif // HISTORY
