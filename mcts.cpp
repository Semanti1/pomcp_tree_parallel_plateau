#include "mcts.h"
#include "testsimulator.h"
#include <math.h>
#include <thread>
#include <algorithm>
#include <atomic>
#include <omp.h>
#include <map>
using namespace std;
using namespace UTILS;

//-----------------------------------------------------------------------------
std::atomic<bool> lock_ = { false };
std::atomic<bool> lock_g = { false };
std::atomic<bool> lock_q = { false };
//void lock() { while (lock_.exchange(true, std::memory_order_acquire)); }
void lockq() { while (lock_q.exchange(true, std::memory_order_acquire)); }
void lockg() { while (lock_g.exchange(true, std::memory_order_acquire)); }
std::map<VNODE*, std::mutex*> vnodelocks;
std::map<QNODE*, std::mutex*> qnodelocks;
//void unlock() { lock_.store(false, std::memory_order_release); }
void unlockq() { lock_q.store(false, std::memory_order_release); }
void unlockg() { lock_g.store(false, std::memory_order_release); }
std::mutex vmtx;
std::mutex qmtx;
std::mutex gmtx;
std::mutex hmtx;

MCTS::PARAMS::PARAMS()
:   Verbose(0),
    MaxDepth(100),
    NumSimulations(1000),
    NumStartStates(1000),
    UseTransforms(true),
    NumTransforms(0),
    MaxAttempts(0),
    ExpandCount(1),
    ExplorationConstant(1),
    UseRave(false),
    RaveDiscount(1.0),
    RaveConstant(0.01),
    DisableTree(false)
{
}

MCTS::MCTS(const SIMULATOR& simulator, const PARAMS& params)
:   Simulator(simulator),
    Params(params),
    TreeDepth(0)
{
    VNODE::NumChildren = Simulator.GetNumActions();
    QNODE::NumChildren = Simulator.GetNumObservations();

    Root = ExpandNode(Simulator.CreateStartState(),History);

    for (int i = 0; i < Params.NumStartStates; i++)
        Root->Beliefs().AddSample(Simulator.CreateStartState());
}

MCTS::~MCTS()
{
    VNODE::Free(Root, Simulator);
    VNODE::FreeAll();
}

bool MCTS::Update(int action, int observation, double reward)
{
    History.Add(action, observation);
    BELIEF_STATE beliefs;

    // Find matching vnode from the rest of the tree
    QNODE& qnode = Root->Child(action);
    VNODE* vnode = qnode.Child(observation);
    if (vnode)
    {
        if (Params.Verbose >= 1)
            cout << "Matched " << vnode->Beliefs().GetNumSamples() << " states" << endl;
        beliefs.Copy(vnode->Beliefs(), Simulator);
    }
    else
    {
        if (Params.Verbose >= 1)
            cout << "No matching node found" << endl;
    }

    // Generate transformed states to avoid particle deprivation
    if (Params.UseTransforms)
        AddTransforms(Root, beliefs);

    // If we still have no particles, fail
    if (beliefs.Empty() && (!vnode || vnode->Beliefs().Empty()))
        return false;

    if (Params.Verbose >= 1)
        Simulator.DisplayBeliefs(beliefs, cout);

    // Find a state to initialise prior (only requires fully observed state)
    const STATE* state = 0;
    if (vnode && !vnode->Beliefs().Empty())
        state = vnode->Beliefs().GetSample(0);
    else
        state = beliefs.GetSample(0);

    // Delete old tree and create new root
    VNODE::Free(Root, Simulator);
    VNODE* newRoot = ExpandNode(state,History);
    newRoot->Beliefs() = beliefs;
    Root = newRoot;
    return true;
}

int MCTS::SelectAction()
{
    if (Params.DisableTree)
        RolloutSearch();
    else
        UCTSearch();
    return GreedyUCB(Root, false);
}

void MCTS::RolloutSearch()
{
	/*std::vector<double> totals(Simulator.GetNumActions(), 0.0);
	int historyDepth = History.Size();
	std::vector<int> legal;
	assert(BeliefState().GetNumSamples() > 0);
	Simulator.GenerateLegal(*BeliefState().GetSample(0), GetHistory(), legal, GetStatus());
	random_shuffle(legal.begin(), legal.end());

	for (int i = 0; i < Params.NumSimulations; i++)
	{
		int action = legal[i % legal.size()];
		STATE* state = Root->Beliefs().CreateSample(Simulator);
		Simulator.Validate(*state);

		int observation;
		double immediateReward, delayedReward, totalReward;
		bool terminal = Simulator.Step(*state, action, observation, immediateReward);

		VNODE*& vnode = Root->Child(action).Child(observation);
		if (!vnode && !terminal)
		{
			vnode = ExpandNode(state);
			AddSample(vnode, *state);
		}
		History.Add(action, observation);

		delayedReward = Rollout(*state,TreeDepth);
		totalReward = immediateReward + Simulator.GetDiscount() * delayedReward;
		Root->Child(action).Value.Add(totalReward);

		Simulator.FreeState(state);
		History.Truncate(historyDepth);
	}*/
}

void MCTS::UCTSearch()
{
    ClearStatistics();
    vnodelocks = {};
    qnodelocks = {};
    //cout << "vnodelocks  "<< vnodelocks.size() << endl;
    //cout << "qnodelocks  " << qnodelocks.size() << endl;
    int historyDepth = History.Size();
    //TreeDepth = 0;
    //PeakTreeDepth = 0;
    omp_set_num_threads(6);
    vnodelocks.emplace(Root, new std::mutex);
    vmtx.lock();
    auto lockithereR = vnodelocks[Root];
    vmtx.unlock();
    #pragma omp parallel for //shared(Root,Simulator,Params,History) //firstprivate(TreeDepth,PeakTreeDepth) //reduction(max:TreeDepth,max:PeakTreeDepth)
    for (int n = 0; n < Params.NumSimulations; n++)
    {
       /* int tid = omp_get_thread_num();
        if (tid == 0) {
            int nthreads = omp_get_num_threads();
            printf("Number of threads = %d\n", nthreads);
        }*/
        STATE* state = 0;
        HISTORY hthread = History;
        

        /*vmtx.lock();
        auto lockithereR = vnodelocks[Root];
        vmtx.unlock();*/
           
        state = Root->Beliefs().CreateSample(Simulator, lockithereR);
        //state = Root->Beliefs().CreateSample(Simulator, Root->getmutex());
            Simulator.Validate(*state);
      
        Status.Phase = SIMULATOR::STATUS::TREE;
        
       
        int treedepth_p = 0;
        int peaktreedepth_p = 0;
        double totalReward = SimulateV(*state, Root, treedepth_p, peaktreedepth_p, hthread);
        
        if (Params.Verbose >= 3)
            DisplayValue(4, cout);

     
        Simulator.FreeState(state);

    }
   

    DisplayStatistics(cout);
}

double MCTS::SimulateV(STATE& state, VNODE* vnode, int& treedepth, int& peaktreedepth,HISTORY& hist)
{
   // lock();
    vmtx.lock();
    auto lockit = vnodelocks.find(vnode);
    if (!(lockit != vnodelocks.end()))
    {
        vnodelocks.emplace(vnode, new std::mutex);
    }
    //vnode->Value.Addm(-1, vnodelocks[vnode]);
    //vnode->Value.Virtualloss(-1);
    vmtx.unlock();
    int action;

    { action = GreedyUCB(vnode, true); }

    //PeakTreeDepth = TreeDepth;
    peaktreedepth = treedepth;
    
    if (treedepth >= Params.MaxDepth) // search horizon reached
        return 0;
   
    if (treedepth == 1)
    {
       
        
            AddSample(vnode, state);
        
      
    }
   
    QNODE& qnode = vnode->Child(action);
    
    qmtx.lock();
    QNODE* newqnode = &qnode;
   // auto lockitq = qnodelocks.find(newqnode);
    auto lockitq = qnodelocks.find(&qnode);
    if (!(lockitq != qnodelocks.end()))
    {
        //qnodelocks.emplace(newqnode, new std::mutex);
        qnodelocks.emplace(&qnode, new std::mutex);
    }
    //qnode.Value.Addm(-1, qnodelocks[&qnode]);
    //qnode.Value.Virtualloss(-1);
    qmtx.unlock();
   
    double totalReward = SimulateQ(state, qnode, action, treedepth,peaktreedepth,hist);
  
    vmtx.lock();
    auto lockithere = vnodelocks[vnode];
    vmtx.unlock();
    vnode->Value.Addm(totalReward, lockithere);
    //vnode->Value.Addm(totalReward);
    //vnode->Value.Addm(totalReward, vnode->getmutex());
    AddRave(vnode, totalReward, treedepth, hist,lockithere);
    //AddRave(vnode, totalReward, treedepth, hist, vnode->getmutex());
   
    return totalReward;
}

double MCTS::SimulateQ(STATE& state, QNODE& qnode, int action, int& treedepth, int& peaktreedepth, HISTORY& hist)
{
    int observation;
    double immediateReward, delayedReward = 0;

    if (Simulator.HasAlpha())
        Simulator.UpdateAlpha(qnode, state);
    bool terminal = Simulator.Step(state, action, observation, immediateReward);
    assert(observation >= 0 && observation < Simulator.GetNumObservations());
//#pragma omp critical
    {
       // History.Add(action, observation);
        hist.Add(action, observation);

    }
  

    if (Params.Verbose >= 3)
    {
        Simulator.DisplayAction(action, cout);
        Simulator.DisplayObservation(state, observation, cout);
        Simulator.DisplayReward(immediateReward, cout);
        Simulator.DisplayState(state, cout);
    }
   // lock();
    int ct = 0;
    VNODE*& vnode = qnode.Child(observation);
    /*vmtx.lock();
    auto lockit = vnodelocks.find(vnode);
    if (!(lockit != vnodelocks.end()))
    {
        vnodelocks.emplace(vnode, new std::mutex);
    }
    //vnode->Value.Addm(-1, vnodelocks[vnode]);
    //vnode->Value.Virtualloss(-1);
    auto lockvnode = vnodelocks[vnode];
    vmtx.unlock();*/
    //auto lockvnode = vnodelocks[vnode];
   // (*lockvnode).lock();
    if (!vnode && !terminal && qnode.Value.GetCount() >= Params.ExpandCount)
        vnode = ExpandNode(&state,hist);

   //(*lockvnode).unlock();

    if (!terminal)
    {
        //#pragma omp atomic
        //TreeDepth++;
        treedepth = treedepth + 1;
        if (vnode)
            delayedReward = SimulateV(state, vnode, treedepth,peaktreedepth,hist);
        else
            delayedReward = Rollout(state, treedepth, hist);
        //#pragma omp atomic
        //TreeDepth--;
        treedepth = treedepth - 1;
    }

    double totalReward = immediateReward + Simulator.GetDiscount() * delayedReward;
    //lock();
    qmtx.lock();
   
    auto lockithereq = qnodelocks[&qnode];
    //qnode.Value.Addm(totalReward, qnodelocks[&qnode]);
    qmtx.unlock();
    qnode.Value.Addm(totalReward, lockithereq);
    //qnode.Value.Addm(totalReward);
    //qnode.Value.Add(totalReward);
   


    return totalReward;
}

void MCTS::AddRave(VNODE* vnode, double totalReward, int treedepth, HISTORY hist, std::mutex* m)
{
   std::lock_guard<std::mutex> lock(*m);
   /*vmtx.lock();        
    auto lockithere = vnodelocks[vnode];
    vmtx.unlock();
    //std::lock_guard<std::mutex> lock(*(vnodelocks[vnode]));
    std::lock_guard<std::mutex> lock(*lockithere);*/
    double totalDiscount = 1.0;
   
    
        for (int t = treedepth; t < hist.Size(); ++t)
    {
       
            QNODE& qnode = vnode->Child(hist[t].Action);
          
           qnode.AMAF.Add(totalReward, totalDiscount);

        totalDiscount *= Params.RaveDiscount;
      //  unlockg();
    }
      
}

VNODE* MCTS::ExpandNode(const STATE* state, HISTORY& hist)
{
    VNODE* vnode = VNODE::Create();
    vnode->Value.Set(0, 0);
   
    Simulator.Prior(state, hist, vnode, Status);

    if (Params.Verbose >= 2)
    {
        cout << "Expanding node: ";
        History.Display(cout);
        cout << endl;
    }

    return vnode;
}

void MCTS::AddSample(VNODE* node, const STATE& state)
{
    
    STATE* sample = Simulator.Copy(state);
    //lock();
    vmtx.lock();
    auto lockithere = vnodelocks[node];
    vmtx.unlock();
    
   // node->Beliefs().AddSamplem(sample, vnodelocks[node]);
    node->Beliefs().AddSamplem(sample, lockithere);
    //node->Beliefs().AddSamplem(sample, node->getmutex());
    //unlock();
    if (Params.Verbose >= 2)
    {
        cout << "Adding sample:" << endl;
        Simulator.DisplayState(*sample, cout);
    }
}

int MCTS::GreedyUCB(VNODE* vnode, bool ucb) const
{
    //static vector<int> besta;
    vector<int> besta;
    besta.clear();
    double bestq = -Infinity;
    //lock();
    int N = vnode->Value.GetCount();
    //unlock();
    double logN = log(N + 1);

    // Locking needed here and in Simulator.GetNumActions()?
    bool hasalpha = Simulator.HasAlpha();

    for (int action = 0; action < Simulator.GetNumActions(); action++)
    {
        double q, alphaq;
        int n, alphan;
        double qc_amaf, qv_amaf;
        //lock();
        QNODE& qnode = vnode->Child(action);
        
        q = qnode.Value.GetValue();
        n = qnode.Value.GetCount();
        qc_amaf = qnode.AMAF.GetCount();
        qv_amaf = qnode.AMAF.GetValue();
       // unlock();
       // if (Params.UseRave && qnode.AMAF.GetCount() > 0)
        if (Params.UseRave && qc_amaf > 0)
        {
            double n2 = qc_amaf;
            double beta = n2 / (n + n2 + Params.RaveConstant * n * n2);
           // q = (1.0 - beta) * q + beta * qnode.AMAF.GetValue();
            q = (1.0 - beta) * q + beta * qv_amaf;
        }

        if (hasalpha && n > 0)
        {
            Simulator.AlphaValue(qnode, alphaq, alphan);
            q = (n * q + alphan * alphaq) / (n + alphan);
            //cout << "N = " << n << ", alphaN = " << alphan << endl;
            //cout << "Q = " << q << ", alphaQ = " << alphaq << endl;
        }

        if (ucb)
            q += FastUCB(N, n, logN);

        if (q >= bestq)
        {
            if (q > bestq)
                besta.clear();
            bestq = q;
            besta.push_back(action);
        }
    }

    assert(!besta.empty());
    return besta[Random(besta.size())];
}

double MCTS::Rollout(STATE& state, int treedepth, HISTORY& hist)
{
    Status.Phase = SIMULATOR::STATUS::ROLLOUT;
    if (Params.Verbose >= 3)
        cout << "Starting rollout" << endl;

    double totalReward = 0.0;
    double discount = 1.0;
    bool terminal = false;
    int numSteps;
    //for (numSteps = 0; numSteps + TreeDepth < Params.MaxDepth && !terminal; ++numSteps)
    for (numSteps = 0; numSteps + treedepth < Params.MaxDepth && !terminal; ++numSteps)
    {
        int observation;
        double reward;
        int action;
//#pragma omp critical
      //  {
        //action = Simulator.SelectRandom(state, History, Status);
        action = Simulator.SelectRandom(state, hist, Status);
        terminal = Simulator.Step(state, action, observation, reward);
      //  #pragma omp critical
        {
       // History.Add(action, observation);
            hist.Add(action, observation);
       }

        if (Params.Verbose >= 4)
        {
            Simulator.DisplayAction(action, cout);
            Simulator.DisplayObservation(state, observation, cout);
            Simulator.DisplayReward(reward, cout);
            Simulator.DisplayState(state, cout);
        }

        totalReward += reward * discount;
        discount *= Simulator.GetDiscount();
    }
    //lock();
#pragma omp critical
    {
    StatRolloutDepth.Add(numSteps); }
   // unlock();
    if (Params.Verbose >= 3)
        cout << "Ending rollout after " << numSteps
            << " steps, with total reward " << totalReward << endl;
    return totalReward;
}

void MCTS::AddTransforms(VNODE* root, BELIEF_STATE& beliefs)
{
    int attempts = 0, added = 0;

    // Local transformations of state that are consistent with history
    while (added < Params.NumTransforms && attempts < Params.MaxAttempts)
    {
        STATE* transform = CreateTransform();
        if (transform)
        {
            beliefs.AddSample(transform);
            added++;
        }
        attempts++;
    }

    if (Params.Verbose >= 1)
    {
        cout << "Created " << added << " local transformations out of "
            << attempts << " attempts" << endl;
    }
}

STATE* MCTS::CreateTransform() const
{
    int stepObs;
    double stepReward;

    STATE* state = Root->Beliefs().CreateSample(Simulator,vnodelocks[Root]);
    Simulator.Step(*state, History.Back().Action, stepObs, stepReward);
    if (Simulator.LocalMove(*state, History, stepObs, Status))
        return state;
    Simulator.FreeState(state);
    return 0;
}

double MCTS::UCB[UCB_N][UCB_n];
bool MCTS::InitialisedFastUCB = true;

void MCTS::InitFastUCB(double exploration)
{
    cout << "Initialising fast UCB table... ";
    for (int N = 0; N < UCB_N; ++N)
        for (int n = 0; n < UCB_n; ++n)
            if (n == 0)
                UCB[N][n] = Infinity;
            else
                UCB[N][n] = exploration * sqrt(log(N + 1) / n);
    cout << "done" << endl;
    InitialisedFastUCB = true;
}

inline double MCTS::FastUCB(int N, int n, double logN) const
{
    if (InitialisedFastUCB && N < UCB_N && n < UCB_n)
        return UCB[N][n];

    if (n == 0)
        return Infinity;
    else
        return Params.ExplorationConstant * sqrt(logN / n);
}

void MCTS::ClearStatistics()
{
    StatTreeDepth.Clear();
    StatRolloutDepth.Clear();
    StatTotalReward.Clear();
}

void MCTS::DisplayStatistics(ostream& ostr) const
{
    if (Params.Verbose >= 1)
    {
        StatTreeDepth.Print("Tree depth", ostr);
        StatRolloutDepth.Print("Rollout depth", ostr);
        StatTotalReward.Print("Total reward", ostr);
    }

    if (Params.Verbose >= 2)
    {
        ostr << "Policy after " << Params.NumSimulations << " simulations" << endl;
        DisplayPolicy(6, ostr);
        ostr << "Values after " << Params.NumSimulations << " simulations" << endl;
        DisplayValue(6, ostr);
    }
}

void MCTS::DisplayValue(int depth, ostream& ostr) const
{
    HISTORY history;
    ostr << "MCTS Values:" << endl;
    Root->DisplayValue(history, depth, ostr);
}

void MCTS::DisplayPolicy(int depth, ostream& ostr) const
{
    HISTORY history;
    ostr << "MCTS Policy:" << endl;
    Root->DisplayPolicy(history, depth, ostr);
}

//-----------------------------------------------------------------------------

void MCTS::UnitTest()
{
    UnitTestGreedy();
    UnitTestUCB();
    UnitTestRollout();
    for (int depth = 1; depth <= 3; ++depth)
        UnitTestSearch(depth);
}

void MCTS::UnitTestGreedy()
{
    TEST_SIMULATOR testSimulator(5, 5, 0);
    PARAMS params;
    MCTS mcts(testSimulator, params);
    int numAct = testSimulator.GetNumActions();
    int numObs = testSimulator.GetNumObservations();

    VNODE* vnode = mcts.ExpandNode(testSimulator.CreateStartState(),mcts.History);
    vnode->Value.Set(1, 0);
    vnode->Child(0).Value.Set(0, 1);
    for (int action = 1; action < numAct; action++)
        vnode->Child(action).Value.Set(0, 0);
    assert(mcts.GreedyUCB(vnode, false) == 0);
}

void MCTS::UnitTestUCB()
{
    TEST_SIMULATOR testSimulator(5, 5, 0);
    PARAMS params;
    MCTS mcts(testSimulator, params);
    int numAct = testSimulator.GetNumActions();
    int numObs = testSimulator.GetNumObservations();

    // With equal value, action with lowest count is selected
    VNODE* vnode1 = mcts.ExpandNode(testSimulator.CreateStartState(),mcts.History);
    vnode1->Value.Set(1, 0);
    for (int action = 0; action < numAct; action++)
        if (action == 3)
            vnode1->Child(action).Value.Set(99, 0);
        else
            vnode1->Child(action).Value.Set(100 + action, 0);
    assert(mcts.GreedyUCB(vnode1, true) == 3);

    // With high counts, action with highest value is selected
    VNODE* vnode2 = mcts.ExpandNode(testSimulator.CreateStartState(), mcts.History);
    vnode2->Value.Set(1, 0);
    for (int action = 0; action < numAct; action++)
        if (action == 3)
            vnode2->Child(action).Value.Set(99 + numObs, 1);
        else
            vnode2->Child(action).Value.Set(100 + numAct - action, 0);
    assert(mcts.GreedyUCB(vnode2, true) == 3);

    // Action with low value and low count beats actions with high counts
    VNODE* vnode3 = mcts.ExpandNode(testSimulator.CreateStartState(),mcts.History);
    vnode3->Value.Set(1, 0);
    for (int action = 0; action < numAct; action++)
        if (action == 3)
            vnode3->Child(action).Value.Set(1, 1);
        else
            vnode3->Child(action).Value.Set(100 + action, 1);
    assert(mcts.GreedyUCB(vnode3, true) == 3);

    // Actions with zero count is always selected
    VNODE* vnode4 = mcts.ExpandNode(testSimulator.CreateStartState(), mcts.History);
    vnode4->Value.Set(1, 0);
    for (int action = 0; action < numAct; action++)
        if (action == 3)
            vnode4->Child(action).Value.Set(0, 0);
        else
            vnode4->Child(action).Value.Set(1, 1);
    assert(mcts.GreedyUCB(vnode4, true) == 3);
}

void MCTS::UnitTestRollout()
{
    TEST_SIMULATOR testSimulator(2, 2, 0);
    PARAMS params;
    params.NumSimulations = 1000;
    params.MaxDepth = 10;
    MCTS mcts(testSimulator, params);
    double totalReward;
    for (int n = 0; n < mcts.Params.NumSimulations; ++n)
    {
        STATE* state = testSimulator.CreateStartState();
        mcts.TreeDepth = 0;
        totalReward += mcts.Rollout(*state, mcts.TreeDepth,mcts.History);
    }
    double rootValue = totalReward / mcts.Params.NumSimulations;
    double meanValue = testSimulator.MeanValue();
    assert(fabs(meanValue - rootValue) < 0.1);
}

void MCTS::UnitTestSearch(int depth)
{
    TEST_SIMULATOR testSimulator(3, 2, depth);
    PARAMS params;
    params.MaxDepth = depth + 1;
    params.NumSimulations = pow(10, depth + 1);
    MCTS mcts(testSimulator, params);
    mcts.UCTSearch();
    double rootValue = mcts.Root->Value.GetValue();
    double optimalValue = testSimulator.OptimalValue();
    assert(fabs(optimalValue - rootValue) < 0.1);
}

//-----------------------------------------------------------------------------
