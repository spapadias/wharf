#ifndef CONFIG_H
#define CONFIG_H

#define dygrl dynamic_graph_representation_learning_with_metropolis_hastings

#include <utility.h>
#include <types.h>
#include <globals.h>

auto graph_update_time_on_insert = timer("GraphUpdateTimeOnInsert", false);
auto walk_update_time_on_insert  = timer("WalkUpdateTimeOnInsert", false);
auto Walking_new_sampling_time   = timer("WalkInsertTimeFor2Jobs", false);
auto Walking_insert_new_samples  = timer("WalkInsertTimeFor2Accs", false);
auto MAV_time                    = timer("MAVTime", false);

auto graph_update_time_on_delete = timer("GraphUpdateTimeOnDelete", false);
auto walk_update_time_on_delete  = timer("WalkUpdateTimeOnDelete", false);

auto merge_calc_triplets_to_delete = timer("MergeCalcTripletsToDelete", false);
auto merge_create_delete_walks     = timer("CreateDeleteWalks", false);
auto merge_multiinsert_ctress      = timer("MergeMultiinsertCtrees", false);
auto Merge_time                    = timer("MergeAllTimer", false);
auto LastMerge                     = timer("LastMerge", false);
auto ReadWalks                     = timer("ReadWalks", false);

// Min and Max Measurements
auto MAV_min = 1000.0;
auto MAV_max = 0.0;
auto Merge_min = 1000.0;
auto Merge_max = 0.0;
auto WalkSampling_min = 1000.0;
auto WalkSampling_max = 0.0;
auto WalkInsert_min = 1000.0;
auto WalkInsert_max = 0.0;

#endif
