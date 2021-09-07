#ifndef NODE_H
#define NODE_H

#include "beliefstate.h"
#include "utils.h"
#include <iostream>
#include<shared_mutex>
#include<mutex>
#include <omp.h>
#include <atomic>
class HISTORY;
class SIMULATOR;
class QNODE;
class VNODE;

//-----------------------------------------------------------------------------
// Efficient computation of value from alpha vectors
// Only used for explicit POMDPs
struct ALPHA
{
    std::vector<double> AlphaSum;
    double MaxValue;
};

//-----------------------------------------------------------------------------

template<class COUNT>
class VALUE
{
public:

    void Set(double count, double value)
    {
      //  std::lock_guard<std::mutex> lock(m);
        //std::unique_lock <std::shared_mutex> lock(m);
//#pragma omp critical
        {    Count = count;
        //#pragma omp atomic
        Total = value * count; }

        
    }

    void Add(double totalReward )
    {
      //  std::lock_guard<std::mutex> lock(m);
        //std::unique_lock <std::shared_mutex> lock(m);
//#pragma omp atomic
        Count += 1.0;
//#pragma omp atomic
        Total += totalReward;
    }

    void Addm(double totalReward, std::mutex* m)
    {
       // m->lock();
        std::lock_guard<std::mutex> lock(*m);
        //std::unique_lock <std::shared_mutex> lock(m);
//#pragma omp atomic
        Count += 1.0;
        //#pragma omp atomic
        Total += totalReward;
     //   m->unlock();
    }
    void Virtualloss(double loss)
    {
       // std::lock_guard<std::mutex> lock(m);
      //  std::unique_lock <std::shared_mutex> lock(m);
//#pragma omp atomic
        Total = Total +  loss;
    }

    void Add(double totalReward, COUNT weight)
    {
       // std::lock_guard<std::mutex> lock(m);
        //std::unique_lock <std::shared_mutex> lock(m);
//#pragma omp atomic
        Count += weight;
//#pragma omp atomic
        Total += totalReward * weight;
    }

    void Addm(double totalReward, COUNT weight, std::mutex* m)
    {
         std::lock_guard<std::mutex> lock(*m);
       //  std::unique_lock <std::shared_mutex> lock(m);
 //#pragma omp atomic
       // m->lock();
        Count += weight;
        //#pragma omp atomic
        Total += totalReward * weight;
       // m->unlock();
    }

    double GetValue() const
    {
      //  std::lock_guard<std::mutex> lock(m);
      //  std::shared_lock <std::shared_mutex> lock(m);
//#pragma omp atomic
        return Count == 0 ? Total : Total / Count;
    }

    COUNT GetCount() const
    {
      //  std::lock_guard<std::mutex> lock(m);
       // std::shared_lock <std::shared_mutex> lock(m);
        return Count;
    }

    /*VALUE operator=(const VALUE& other)
    {
        if (this != &other) {
            std::unique_lock<std::shared_mutex> _mylock(m, std::defer_lock);
            std::shared_lock<std::shared_mutex> _otherlock(other.m, std::defer_lock);
            std::lock(_mylock, _otherlock);
            Total = other.Total;
            Count = other.Count;
        }
        return *this;
    }*/

private:

    COUNT Count;
    double Total;
    //std::atomic<COUNT> Count;
   // std::atomic<double> Total;
   // mutable  std::shared_mutex m;
};

//-----------------------------------------------------------------------------

class QNODE
{
public:

    VALUE<int> Value;
    VALUE<double> AMAF;

    void Initialise();

    VNODE*& Child(int c) { return Children[c]; }
    VNODE* Child(int c) const { return Children[c]; }
    ALPHA& Alpha() { return AlphaData; }
    const ALPHA& Alpha() const { return AlphaData; }

    void DisplayValue(HISTORY& history, int maxDepth, std::ostream& ostr) const;
    void DisplayPolicy(HISTORY& history, int maxDepth, std::ostream& ostr) const;
    static int getnumchildren() { return NumChildren; }
    static int NumChildren;
    //int NumChildren;

private:

    std::vector<VNODE*> Children;
    ALPHA AlphaData;
  //  std::mutex m;

friend class VNODE;
};

//-----------------------------------------------------------------------------

class VNODE : public MEMORY_OBJECT
{
public:

    VALUE<int> Value;

    void Initialise();
    static VNODE* Create();
   // VNODE* Create();
    static void Free(VNODE* vnode, const SIMULATOR& simulator);
    static void FreeAll();

    QNODE& Child(int c) {// std::lock_guard<std::mutex> lock(mchild); 
        return Children[c]; }
    const QNODE& Child(int c) const {  return Children[c]; }
    BELIEF_STATE& Beliefs() {  return BeliefState; }
    const BELIEF_STATE& Beliefs() const {  return BeliefState; }
    //std::mutex* getmutex() { return mchild; }
    void SetChildren(int count, double value);

    void DisplayValue(HISTORY& history, int maxDepth, std::ostream& ostr) const;
    void DisplayPolicy(HISTORY& history, int maxDepth, std::ostream& ostr) const;

    static int NumChildren;
    //int NumChildren;

private:

    std::vector<QNODE> Children;
   // std::unique_ptr<std::vector<QNODE>> Children;
    BELIEF_STATE BeliefState;
    static MEMORY_POOL<VNODE> VNodePool;
    //MEMORY_POOL<VNODE> VNodePool;
  // std::mutex* mchild;
  //  std::mutex mbeliefst;
  //  std::mutex mvnodepool;
};

#endif // NODE_H
