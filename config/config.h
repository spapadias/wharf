#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_CONFIG_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_CONFIG_H

#define dygrl dynamic_graph_representation_learning_with_metropolis_hastings

#include <utility.h>
#include <types.h>
#include <globals.h>

// walk and graph update times
auto graph_update_timer_insert = timer("GraphUpdateTimerInsert", false);
auto walk_update_timer_insert  = timer("WalkUpdateTimerInsert", false);

auto graph_update_timer_delete = timer("GraphUpdateTimerDelete", false);
auto walk_update_timer_delete  = timer("WalkUpdateTimerDelete", false);

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_CONFIG_H
