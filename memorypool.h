#ifndef MEMORY_POOL_H
#define MEMORY_POOL_H

#include <vector>
#include <ostream>
#include<mutex>
class MEMORY_OBJECT
{
public:

    void SetAllocated() { //std::lock_guard<std::mutex> lock(mutex);
        Allocated = true; }
    void ClearAllocated() {// std::lock_guard<std::mutex> lock(mutex);
        Allocated = false; }
    bool IsAllocated() const {  return Allocated; }

private:

    bool Allocated;
   // mutable std::mutex mutex;
};

template <class T>
class MEMORY_POOL
{
public:

    MEMORY_POOL()
    : NumAllocated(0)
    {
    }

    ~MEMORY_POOL()
    {
        DeleteAll();
    }

    T* Construct()
    {
        T* obj = Allocate();
        return new (obj) T;
    }

    void Destroy(T* obj)
    {
        obj.T::~T();
        Free(obj);
    }

    T* Allocate() 
    { 
        std::lock_guard<std::mutex> lock(mutex_freelist);
        //std::scoped_lock<std::mutex, std::mutex> lock(mutex_freelist, mutex_numallocated);
        /*std::lock(m1, m2);
        std::lock_guard<std::mutex> lck1(m1, std::adopt_lock);
        std::lock_guard<std::mutex> lck2(m2, std::adopt_lock);*/
        if (FreeList.empty())
            NewChunk();
        T* obj = FreeList.back();
        FreeList.pop_back();
        assert(!obj->IsAllocated());
        obj->SetAllocated();
        NumAllocated++;
        return obj;
    }
    
    void Free(T* obj)
    { 
       std::lock_guard<std::mutex> lock(mutex_freelist);
       // std::scoped_lock<std::mutex, std::mutex> lock(mutex_freelist, mutex_numallocated);
        assert(obj->IsAllocated());
        obj->ClearAllocated();
        FreeList.push_back(obj);
        NumAllocated--;
    }
    
    void DeleteAll()
    {
        std::lock_guard<std::mutex> lock(mutex_freelist);
        //std::scoped_lock<std::mutex, std::mutex, std::mutex> lock(mutex_freelist, mutex_numallocated, mutex_chunks);
        for (ChunkIterator i_chunk = Chunks.begin(); i_chunk != Chunks.end(); ++i_chunk)
            delete *i_chunk;
        Chunks.clear();
        FreeList.clear();
        NumAllocated = 0;
    }
    
    int GetNumAllocated() const { return NumAllocated; }

private:

    struct CHUNK
    {
        static const int Size = 256;
        T Objects[Size];
    };

    void NewChunk()
    {
        //std::scoped_lock<std::mutex, std::mutex> lock(mutex_freelist, mutex_chunks);
        //std::lock_guard<std::mutex> lock(mutex_freelist);
        CHUNK* chunk = new CHUNK;
        Chunks.push_back(chunk);
        for (int i = CHUNK::Size - 1; i >= 0; --i)
        {
            FreeList.push_back(&chunk->Objects[i]);
            chunk->Objects[i].ClearAllocated();
        }
    }

    std::vector<CHUNK*> Chunks;
    std::mutex mutex_freelist;
    //std::mutex mutex_chunks;
    std::vector<T*> FreeList;
    int NumAllocated;
   // std::mutex mutex_numallocated;
    typedef typename std::vector<CHUNK*>::iterator ChunkIterator;
};

#endif // MEMORY_POOL_H
