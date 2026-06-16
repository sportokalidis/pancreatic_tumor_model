/*
 * fake_numa.c — LD_PRELOAD shim for BioDynaMo in Singularity containers.
 *
 * Problem: this BioDynaMo build requires NUMA (fatals if numa_available()==-1)
 * but numa_alloc_onnode() fails on node 1 inside containers with cgroup memory
 * restrictions (multi-NUMA HPC nodes).
 *
 * Fix: report NUMA as available (return 0) but advertise only 1 configured
 * node. BioDynaMo then assigns all threads to domain 0 and never calls
 * RunOnNuma() for domain 1, so the failing allocation never happens.
 *
 * Build (once, on login node):
 *   gcc -shared -fPIC -nostartfiles -o scripts/hpc/fake_numa.so \
 *       scripts/hpc/fake_numa.c
 */
#include <stddef.h>
#include <stdlib.h>

/* NUMA is "available" — required or BioDynaMo fatals */
int numa_available(void) { return 0; }

/* Advertise 1 NUMA node → BioDynaMo never tries node 1 */
int numa_num_configured_nodes(void) { return 1; }
int numa_max_node(void) { return 0; }
int numa_num_task_cpus(void) { return 1; }
int numa_num_task_nodes(void) { return 1; }

/* All CPUs belong to node 0 — prevents NumaPoolAllocator from indexing
 * into a pool slot that was never allocated (pool is sized for 1 node). */
int numa_node_of_cpu(int cpu) { return 0; }

/* Thread / memory policy stubs */
int  numa_run_on_node(int node)          { return 0; }
void numa_set_localalloc(void)           { }
void numa_set_strict(int flag)           { }
void numa_set_preferred(int node)        { }

/* Redirect all node-specific allocations to plain malloc */
void *numa_alloc_onnode(size_t size, int node) { return malloc(size); }
void *numa_alloc_local(size_t size)             { return malloc(size); }
void *numa_alloc(size_t size)                   { return malloc(size); }
void  numa_free(void *start, size_t size)       { free(start); }
