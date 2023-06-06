#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "path_with_turn_forbidden.h"
#include "floyd_warshall.h"
#include "omp.h"
#include "utils.h"

PathNode *new_path_node(int id, int depth, PathNode *father, int children_num)
{
    PathNode *node = (PathNode *)calloc(1, sizeof(PathNode));
    assert(node != NULL);
    node->id = id;
    node->depth = depth;
    node->father = father;
    node->children_num = children_num;
    return node;
}

void delete_path_line(PathNode *node)
{
    PathNode *father = NULL;
    while (node->children_num == 0)
    {
        father = node->father;
        free(node);
        if (father != NULL)
        {
            father->children_num--;
            node = father;
        }
        else
            break;
    }
}

// shortest path with turn forbidden
void trun_model_dijkstra(const GraphAdjList *G, const TreeCoordinate *tc, int root, unsigned int *dist, int **path, int (*turn_check)(int, int, int, int, int, int))
{
    assert(dist != NULL && path != NULL && tc != NULL && dist != NULL);
    memset(dist, INF, G->num_nodes * sizeof(int));
    int **visit = (int **)calloc(G->num_nodes, sizeof(int *)); // wheather one path is visited
    assert(visit != NULL);
    for (int i = 0; i < G->num_nodes; ++i)
    {
        visit[i] = (int *)calloc(G->num_nodes, sizeof(int));
        assert(visit[i] != NULL);
        memset(path[i], -1, G->num_nodes * sizeof(int));
    }
    PathNode **queue = (PathNode **)calloc(G->num_nodes * G->degree, sizeof(PathNode *));
    int queue_s = 0; // strat
    int queue_e = 0; // end
    // push root in queue
    queue[queue_e++] = new_path_node(root, 0, NULL, 0); // root node
    dist[root] = 0;
    path[root][root] = root;
    // start BFS
    PathNode *this_node = 0;
    PathNode *pre_node = 0;
    // PathNode *dst_node = 0;
    int this_node_id = 0;
    int pre_node_id = 0;
    int dst_node_id = 0;
    int end = 0;
    while (queue_e > queue_s)
    {
        end = queue_e;
        while (queue_s < end)
        {
            this_node = queue[queue_s++];
            this_node_id = this_node->id;
            if (this_node->father != NULL)
            {
                pre_node = this_node->father;
                pre_node_id = pre_node->id;
            }
            else
                pre_node_id = this_node_id;
            for (int i = 0; i < G->degree; ++i)
            {
                dst_node_id = G->adj_list[this_node_id][i];
                if (dst_node_id != INF && visit[this_node_id][dst_node_id] == 0) // not visited path
                {
                    if ((*turn_check)(tc->x[pre_node_id], tc->y[pre_node_id], tc->x[this_node_id], tc->y[this_node_id], tc->x[dst_node_id], tc->y[dst_node_id]) != 0)
                    {
                        visit[this_node_id][dst_node_id] = 1;
                        queue[queue_e++] = new_path_node(dst_node_id, this_node->depth + 1, this_node, 0);
                        assert(queue_e < G->num_nodes * G->degree);
                        this_node->children_num++;
                        if (dist[dst_node_id] == INF)
                        {
                            dist[dst_node_id] = this_node->depth + 1;
                            // record path
                            path[dst_node_id][dst_node_id] = this_node_id;
                            PathNode *leaf_node = this_node;
                            while (leaf_node->father != NULL)
                            {
                                path[dst_node_id][leaf_node->id] = leaf_node->father->id;
                                leaf_node = leaf_node->father;
                            }
                        }
                    }
                }
            }
            if (this_node->children_num == 0)
                delete_path_line(this_node);
        }
    }
    free(queue);
    for (int i = 0; i < G->num_nodes; ++i)
        free(visit[i]);
    free(visit);
}

void trun_model_dijkstra_parallel(const GraphAdjList *G, const TreeCoordinate *tc, unsigned int **dist, int ***path, int (*turn_check)(int, int, int, int, int, int))
{
    omp_set_num_threads(PROCS_NUM);
#pragma omp parallel for
    for (int i = 0; i < G->num_nodes; ++i)
        trun_model_dijkstra(G, tc, i, dist[i], path[i], (*turn_check));
    floyd_check_dist(G, dist);
}

void turn_model_check_path(int num_nodes, int ***path, TreeCoordinate *tc, int (*turn_check)(int, int, int, int, int, int))
{
    // printf("------------------ turn check ---------------------\n");
    omp_set_num_threads(PROCS_NUM);
#pragma omp parallel for
    // turn model check
    for (int i = 0; i < num_nodes; ++i) // src node
    {
        // printf("******* src %d  ********\n", i);
        // print_matrix(path[i], num_nodes, "123");
        int pre_node, this_node, post_node;
        for (int j = 0; j < num_nodes; ++j) // dst node
        {
            post_node = j;
            this_node = path[i][j][post_node];
            while (this_node != i)
            {
                pre_node = path[i][j][this_node];
                assert((*turn_check)(tc->x[pre_node], tc->y[pre_node], tc->x[this_node], tc->y[this_node], tc->x[post_node], tc->y[post_node]) != 0);
                post_node = this_node;
                this_node = pre_node;
            }
        }
    }
}

void turn_model_initial_path_dist(const GraphAdjList *G, unsigned int ***dist, int ****path)
{
    (*dist) = (unsigned int **)calloc(G->num_nodes, sizeof(unsigned int *));
    (*path) = (int ***)calloc(G->num_nodes, sizeof(int **));
    assert(*dist != NULL && *path != NULL);
    for (int i = 0; i < G->num_nodes; ++i)
    {
        (*dist)[i] = (unsigned int *)calloc(G->num_nodes, sizeof(unsigned int));
        (*path)[i] = (int **)calloc(G->num_nodes, sizeof(int *));
        assert((*dist)[i] != NULL && (*path)[i] != NULL);
        for (int j = 0; j < G->num_nodes; ++j)
        {
            (*path)[i][j] = (int *)calloc(G->num_nodes, sizeof(int));
            assert((*path)[i][j] != NULL);
        }
    }
}

void turn_model_delete_path_dist(const GraphAdjList *G, unsigned int **dist, int ***path)
{
    for (int i = 0; i < G->num_nodes; ++i)
    {
        for (int j = 0; j < G->num_nodes; ++j)
            free(path[i][j]);
        free(path[i]);
        free(dist[i]);
    }
    free(dist);
    free(path);
}

double turn_model_cost(const GraphAdjList *G, unsigned int **dist, int ***path, unsigned int *history_path_use_array, int update_en)
{
    // average dist
    unsigned long sum = 0;
    for (int i = 0; i < G->num_nodes; ++i)
        for (int j = 0; j < G->num_nodes; ++j)
            sum += dist[i][j];
    double average = (double)((long double)sum / (long double)(G->num_nodes * G->num_nodes));
    // traffic balance
    int **path_use = (int **)calloc(G->num_nodes, sizeof(int *));
    assert(path_use != NULL);
    for (int i = 0; i < G->num_nodes; ++i)
    {
        path_use[i] = (int *)calloc(G->num_nodes, sizeof(int));
        assert(path_use[i] != NULL);
    }
    int pre_node, post_node;
    for (int i = 0; i < G->num_nodes; ++i)
    {
        for (int j = 0; j < G->num_nodes; ++j)
        {
            post_node = j;
            pre_node = path[i][j][post_node];
            while (pre_node != i)
            {
                path_use[pre_node][post_node] += 1;
                post_node = pre_node;
                pre_node = path[i][j][post_node];
            }
        }
    }
    unsigned int *path_use_array = NULL;
    if (update_en == 0)
    {
        path_use_array = (unsigned int *)calloc(G->num_nodes * G->degree, sizeof(unsigned int));
        memcpy(path_use_array, history_path_use_array, G->num_nodes * G->degree * sizeof(unsigned int));
    }
    else if (update_en == 1)
        path_use_array = history_path_use_array;
    else
        assert(0);
    int ptr = 0;
    for (int i = 0; i < G->num_nodes; ++i)
    {
        for (int j = 0; j < G->num_nodes; ++j)
        {
            if (path_use[i][j] > 0)
            {
                assert(ptr < G->num_nodes * G->degree);
                path_use_array[ptr++] += path_use[i][j];
            }
        }
    }
    double var = variance(G->num_nodes * G->degree, path_use_array);
    // free
    for (int i = 0; i < G->num_nodes; ++i)
        free(path_use[i]);
    free(path_use);
    if (update_en == 0)
        free(path_use_array);
    return 1000.0 * average + 1.0 * var;
}

int *root_selection(const GraphAdjList *G, int (*turn_check)(int, int, int, int, int, int), int root_num)
{
    TreeCoordinate *tc;
    int *best_root = (int *)calloc(root_num, sizeof(int));
    assert(best_root != NULL);
    unsigned int **dist = NULL;
    int ***path = NULL;
    turn_model_initial_path_dist(G, &dist, &path);
    unsigned int *history_path_use_array = (unsigned int *)calloc(G->num_nodes * G->degree, sizeof(unsigned int));
    for (int k = 0; k < root_num; ++k) // every root
    {
        double smallest_cost = 0;
        double this_cost = 0;
        for (int i = 0; i < G->num_nodes; ++i)
        {
            tc = gen_tree_coordinate(G, i, 1);
            trun_model_dijkstra_parallel(G, tc, dist, path, turn_check);
            // cost
            this_cost = turn_model_cost(G, dist, path, history_path_use_array, 0);
            if (this_cost < smallest_cost || i == 0)
            {
                smallest_cost = this_cost;
                best_root[k] = i;
            }
            delete_tc(tc);
        }
        // record path-use info
        tc = gen_tree_coordinate(G, best_root[k], 1);
        trun_model_dijkstra_parallel(G, tc, dist, path, turn_check);
        this_cost = turn_model_cost(G, dist, path, history_path_use_array, 1);
        delete_tc(tc);
    }
    turn_model_delete_path_dist(G, dist, path);
    free(history_path_use_array);
    return best_root;
}

// shortest path with turn forbidden
void trun_model_bfs(const GraphAdjList *G, const TreeCoordinate *tc, int root, unsigned int *dist, int **path, int (*turn_check)(int, int, int, int, int, int))
{
    assert(dist != NULL && path != NULL && tc != NULL && dist != NULL);
    memset(dist, INF, G->num_nodes * sizeof(int));
    int *visit = (int *)calloc(G->num_nodes, sizeof(int)); // wheather one node is visited ( all connected nodes are visited )
    assert(visit != NULL);
    for (int i = 0; i < G->num_nodes; ++i)
        memset(path[i], -1, G->num_nodes * sizeof(int));
    PathNode **queue = (PathNode **)calloc(G->num_nodes * G->degree, sizeof(PathNode *));
    int queue_s = 0; // strat
    int queue_e = 0; // end
    int length = 0;
    int max_length = G->num_nodes * G->degree;
    // push root in queue
    queue[queue_e++] = new_path_node(root, 0, NULL, 0); // root node
    length++;
    dist[root] = 0;
    path[root][root] = root;
    // start BFS
    PathNode *this_node = NULL;
    PathNode *pre_node = NULL;
    // PathNode *dst_node = 0;
    int this_node_id = 0;
    int pre_node_id = 0;
    int dst_node_id = 0;
    int end = 0;
    int all_neighbor_dist = 0;
    while (queue_e != queue_s)
    {
        end = queue_e;
        while (queue_s != end)
        {
            this_node = queue[queue_s++];
            queue_s = queue_s % max_length;
            length--;
            this_node_id = this_node->id;
            if (this_node->father != NULL)
            {
                pre_node = this_node->father;
                pre_node_id = pre_node->id;
            }
            else
                pre_node_id = this_node_id;
            for (int i = 0; i < G->degree; ++i)
            {
                dst_node_id = G->adj_list[this_node_id][i];
                if (dst_node_id != INF && visit[dst_node_id] == 0) // not visited path
                {
                    if ((*turn_check)(tc->x[pre_node_id], tc->y[pre_node_id], tc->x[this_node_id], tc->y[this_node_id], tc->x[dst_node_id], tc->y[dst_node_id]) != 0)
                    {
                        // visit[dst_node_id] = 1;
                        queue[queue_e++] = new_path_node(dst_node_id, this_node->depth + 1, this_node, 0);
                        queue_e = queue_e % max_length;
                        length++;
                        assert(length < max_length);
                        this_node->children_num++;
                        if (dist[dst_node_id] == INF)
                        {
                            dist[dst_node_id] = this_node->depth + 1;
                            // record path
                            path[dst_node_id][dst_node_id] = this_node_id;
                            PathNode *leaf_node = this_node;
                            while (leaf_node->father != NULL)
                            {
                                path[dst_node_id][leaf_node->id] = leaf_node->father->id;
                                leaf_node = leaf_node->father;
                            }
                        }
                    }
                }
            }
            all_neighbor_dist = 0;
            for (int i = 0; i < G->degree; ++i)
            {
                if (G->adj_list[this_node_id][i] != INF && dist[G->adj_list[this_node_id][i]] == INF)
                    all_neighbor_dist = 1;
            }
            if (all_neighbor_dist == 0)
                visit[this_node_id] = 1;
            if (this_node->children_num == 0)
                delete_path_line(this_node);
        }
    }
    free(queue);
    free(visit);
}

void trun_model_bfs_parallel(const GraphAdjList *G, const TreeCoordinate *tc, unsigned int **dist, int ***path, int (*turn_check)(int, int, int, int, int, int))
{
    omp_set_num_threads(PROCS_NUM);
#pragma omp parallel for
    for (int i = 0; i < G->num_nodes; ++i)
        trun_model_bfs(G, tc, i, dist[i], path[i], (*turn_check));
    floyd_check_dist(G, dist);
}

void turn_model_check_dist_new(const GraphAdjList *G, unsigned int ***dist)
{
    int false_path_num = 0;
    for (int d = 0; d < G->degree; ++d)
        for (int i = 0; i < G->num_nodes; ++i)
            for (int j = 0; j < G->num_nodes; ++j)
                if (dist[i][d][j] == INF)
                    false_path_num++;
    for (int d = G->degree; d < G->degree + 1; ++d)
        for (int i = 0; i < G->num_nodes; ++i)
            for (int j = 0; j < G->num_nodes; ++j)
                assert(dist[i][d][j] < INF);
    double false_dist_rate = (double)false_path_num / (double)(G->degree * G->num_nodes * G->num_nodes);
    printf("False dist rate is %lf\n", false_dist_rate);
}

// shortest path with turn forbidden (find every shortest path for every input channel)
void trun_model_bfs_new(const GraphAdjList *G, TreeCoordinate **tc, RouterList *R, int (*turn_check)(int, int, int, int, int, int))
{
    Router *router_list = R->rl;
    // initial LUT and SELF_DST
    for (int i = 0; i < R->num_nodes; ++i)
    {
        for (int j = 0; j < ROUTE_SIM_ANYTOPO_VC; ++j)
        {
            router_list[i].lut[j] = (unsigned char **)calloc(R->degree + 1, sizeof(unsigned char *));
            router_list[i].self_dst[j] = (unsigned int **)calloc(R->degree + 1, sizeof(unsigned int *));
            assert(router_list[i].lut[j] != NULL && router_list[i].self_dst[j] != NULL);
            for (int k = 0; k < R->degree + 1; ++k)
            {
                router_list[i].lut[j][k] = (unsigned char *)calloc(R->num_nodes, sizeof(unsigned char));
                router_list[i].self_dst[j][k] = (unsigned int *)calloc(R->num_nodes, sizeof(unsigned int));
                assert(router_list[i].lut[j][k] != NULL && router_list[i].self_dst[j][k] != NULL);
                memset(router_list[i].lut[j][k], 255, R->num_nodes * sizeof(unsigned char));
                memset(router_list[i].self_dst[j][k], INF, R->num_nodes * sizeof(unsigned int));
            }
        }
    }
    // initial LUT and SELF_DST
    for (int vc_idx = 0; vc_idx < ROUTE_SIM_ANYTOPO_VC; ++vc_idx) // for every VC
    {
        for (int root = 0; root < R->num_nodes; ++root) // for every root node
        {
            for (int input_idx = 0; input_idx < G->degree + 1; ++input_idx) // for every input index
            {
                router_list[root].self_dst[vc_idx][input_idx][root] = 0;
                router_list[root].lut[vc_idx][input_idx][root] = 254;
            }
        }
    }
    // initial dst
    unsigned char **dst = (unsigned char **)calloc(R->num_nodes, sizeof(unsigned char *));
    assert(dst != NULL);
    for (int i = 0; i < R->num_nodes; ++i)
    {
        dst[i] = (unsigned char *)calloc(R->num_nodes, sizeof(unsigned char));
        assert(dst[i] != NULL);
        memset(dst[i], 255, R->num_nodes * sizeof(unsigned char));
        for (unsigned char j = 0; j < G->degree; ++j)
        {
            if (G->adj_list[i][j] != INF)
            {
                dst[i][G->adj_list[i][j]] = j; // dst[i][p] = j; means: node i's j-th port is connected to node p;
            }
        }
    }
    //
    int **visit = (int **)calloc(G->num_nodes, sizeof(int *)); // wheather one path is visited
    assert(visit != NULL);
    for (int i = 0; i < G->num_nodes; ++i)
    {
        visit[i] = (int *)calloc(G->num_nodes, sizeof(int));
        assert(visit[i] != NULL);
    }
    PathNode **queue = (PathNode **)calloc(G->num_nodes * G->degree, sizeof(PathNode *));
    int queue_s = 0;
    int queue_e = 0;
    for (int vc_idx = 0; vc_idx < ROUTE_SIM_ANYTOPO_VC; ++vc_idx) // for every VC
    {
        for (int root = 0; root < R->num_nodes; ++root) // for every root node
        {
            for (int dst_idx = 0; dst_idx < R->num_nodes; ++dst_idx) // every dstenation node
            {
                if (dst_idx == root)
                    continue;
                for (int input_idx = 0; input_idx < G->degree + 1; ++input_idx) // for every input index
                {
                    for (int i = 0; i < G->num_nodes; ++i)
                        memset(visit[i], 0, G->num_nodes * sizeof(int));
                    queue_s = 0; // strat
                    queue_e = 0; // end
                    // push root in queue
                    queue[queue_e++] = new_path_node(root, 0, NULL, 0); // root node
                    // set initial visited path
                    if (input_idx < G->degree)
                    {
                        int virtual_src = G->adj_list[root][input_idx];
                        if (virtual_src != INF)
                            visit[virtual_src][root] = 1;
                    }
                    // start BFS
                    PathNode *this_node = 0;
                    PathNode *pre_node = 0;
                    PathNode *back_node_1 = 0;
                    PathNode *back_node_2 = 0;
                    PathNode *back_node_3 = 0;
                    int this_node_id = 0;
                    int pre_node_id = 0;
                    int dst_node_id = 0;
                    int flag_done = 0;
                    while (flag_done == 0 && queue_e > queue_s)
                    {
                        this_node = queue[queue_s++];
                        this_node_id = this_node->id;
                        int this_node_in_port = 0;
                        if (this_node->father != NULL)
                        {
                            pre_node = this_node->father;
                            pre_node_id = pre_node->id;
                            this_node_in_port = dst[this_node_id][pre_node_id];
                        }
                        else // this node is root
                        {
                            assert(this_node_id == root);
                            if (input_idx == G->degree || G->adj_list[this_node_id][input_idx] == INF)
                                pre_node_id = this_node_id; // send from itself
                            else
                                pre_node_id = G->adj_list[this_node_id][input_idx]; // send from input_idx
                            this_node_in_port = input_idx;
                        }
                        if (router_list[this_node_id].self_dst[vc_idx][this_node_in_port][dst_idx] != INF) // the shortest path is found
                        {
                            assert(router_list[this_node_id].lut[vc_idx][this_node_in_port][dst_idx] != 255);
                            if (this_node_id == root)
                                continue;
                            // set back up the whole path's LUT
                            back_node_1 = this_node;
                            back_node_2 = this_node->father;
                            back_node_3 = back_node_2->father;
                            int back_node_2_in_port = 0;
                            int back_node_1_in_port = 0;
                            while (back_node_2 != NULL)
                            {
                                if (back_node_3 == NULL) // set root
                                {
                                    back_node_2_in_port = input_idx;
                                    back_node_1_in_port = dst[back_node_1->id][back_node_2->id];
                                    assert(router_list[root].self_dst[vc_idx][back_node_2_in_port][dst_idx] == INF);
                                    assert(router_list[root].lut[vc_idx][back_node_2_in_port][dst_idx] == 255);
                                    router_list[root].self_dst[vc_idx][back_node_2_in_port][dst_idx] = router_list[back_node_1->id].self_dst[vc_idx][back_node_1_in_port][dst_idx] + back_node_1->depth; // self_dst
                                    router_list[root].lut[vc_idx][back_node_2_in_port][dst_idx] = dst[back_node_2->id][back_node_1->id];                                                                 // LUT
                                    flag_done = 1;
                                    break;
                                }
                                else
                                {
                                    back_node_2_in_port = dst[back_node_2->id][back_node_3->id];
                                    back_node_1_in_port = dst[back_node_1->id][back_node_2->id];
                                    if (router_list[back_node_2->id].self_dst[vc_idx][back_node_2_in_port][dst_idx] == INF)
                                        router_list[back_node_2->id].self_dst[vc_idx][back_node_2_in_port][dst_idx] = router_list[back_node_1->id].self_dst[vc_idx][back_node_1_in_port][dst_idx] + 1; // self_dst
                                    else
                                        assert(router_list[back_node_2->id].self_dst[vc_idx][back_node_2_in_port][dst_idx] == router_list[back_node_1->id].self_dst[vc_idx][back_node_1_in_port][dst_idx] + 1);
                                    if (router_list[back_node_2->id].lut[vc_idx][back_node_2_in_port][dst_idx] == 255)
                                        router_list[back_node_2->id].lut[vc_idx][back_node_2_in_port][dst_idx] = dst[back_node_2->id][back_node_1->id]; // LUT
                                    else
                                        assert(router_list[back_node_2->id].lut[vc_idx][back_node_2_in_port][dst_idx] == dst[back_node_2->id][back_node_1->id]);
                                    back_node_3 = back_node_3->father;
                                    back_node_2 = back_node_2->father;
                                    back_node_1 = back_node_1->father;
                                }
                            }
                        }
                        else // continue BFS to find shortest path
                        {
                            assert(router_list[this_node_id].lut[vc_idx][this_node_in_port][dst_idx] == 255);
                            int shortest_son_id = -1; // the son, who has the shortest distance to dst_idx
                            unsigned int shortest_distance = INF;
                            int all_son_is_know = 1;
                            // whether this node's children's self_dst is not INF
                            for (int i = 0; i < G->degree; ++i)
                            {
                                dst_node_id = G->adj_list[this_node_id][i];
                                if (dst_node_id != INF && visit[this_node_id][dst_node_id] == 0) // it is a not visited path
                                {
                                    if ((*turn_check)(tc[vc_idx]->x[pre_node_id], tc[vc_idx]->y[pre_node_id], tc[vc_idx]->x[this_node_id], tc[vc_idx]->y[this_node_id], tc[vc_idx]->x[dst_node_id], tc[vc_idx]->y[dst_node_id]) != 0)
                                    {
                                        if (router_list[dst_node_id].self_dst[vc_idx][dst[dst_node_id][this_node_id]][dst_idx] < shortest_distance)
                                        {
                                            shortest_distance = router_list[dst_node_id].self_dst[vc_idx][dst[dst_node_id][this_node_id]][dst_idx];
                                            shortest_son_id = dst_node_id;
                                        }
                                        if (router_list[dst_node_id].self_dst[vc_idx][dst[dst_node_id][this_node_id]][dst_idx] == INF)
                                            all_son_is_know = 0;
                                    }
                                }
                            }
                            if (shortest_son_id != -1 && all_son_is_know == 1)
                            {
                                visit[this_node_id][shortest_son_id] = 1;
                                PathNode *post_node = new_path_node(shortest_son_id, this_node->depth + 1, this_node, 0);
                                queue[queue_e++] = post_node;
                                assert(queue_e < G->num_nodes * G->degree);
                                this_node->children_num++;
                            }
                            else
                            {
                                for (int i = 0; i < G->degree; ++i)
                                {
                                    dst_node_id = G->adj_list[this_node_id][i];
                                    if (dst_node_id != INF && visit[this_node_id][dst_node_id] == 0) // it is a not visited path
                                    {
                                        if ((*turn_check)(tc[vc_idx]->x[pre_node_id], tc[vc_idx]->y[pre_node_id], tc[vc_idx]->x[this_node_id], tc[vc_idx]->y[this_node_id], tc[vc_idx]->x[dst_node_id], tc[vc_idx]->y[dst_node_id]) != 0)
                                        {
                                            visit[this_node_id][dst_node_id] = 1;
                                            PathNode *post_node = new_path_node(dst_node_id, this_node->depth + 1, this_node, 0);
                                            queue[queue_e++] = post_node;
                                            assert(queue_e < G->num_nodes * G->degree);
                                            this_node->children_num++;
                                        }
                                    }
                                }
                            }
                        }
                        if (this_node->children_num == 0)
                            delete_path_line(this_node);
                    }
                }
            }
        }
    }
    free(queue);
    for (int i = 0; i < G->num_nodes; ++i)
    {
        free(visit[i]);
        free(dst[i]);
    }
    free(visit);
    free(dst);
}
