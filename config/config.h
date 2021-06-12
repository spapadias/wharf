#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_CONFIG_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_CONFIG_H

#define dygrl dynamic_graph_representation_learning_with_metropolis_hastings

#include <utility.h>
#include <types.h>
#include <globals.h>

auto graph_update_time_on_insert = timer("GraphUpdateTimeOnInsert", false);
auto walk_update_time_on_insert  = timer("WalkUpdateTimeOnInsert", false);

auto graph_update_time_on_delete = timer("GraphUpdateTimeOnDelete", false);
auto walk_update_time_on_delete  = timer("WalkUpdateTimeOnDelete", false);

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_CONFIG_H
