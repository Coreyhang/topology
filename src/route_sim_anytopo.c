#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include <string.h>
#include "route_sim_anytopo.h"
#include "utils.h"
#include "floyd_warshall.h"
#include "dfs.h"
#include "tree_turn.h"
#include "l_turn.h"
#include "octo_turn.h"
#include "topo_gen.h"
#include "path_with_turn_forbidden.h"

RouterList *construct_router_list_from_graphadjlist(GraphAdjList *G)
{
    RouterList *R = (RouterList *)calloc(1, sizeof(RouterList));
    assert(R != NULL);
    R->num_nodes = G->num_nodes;
    R->degree = G->degree;
    R->rl = (Router *)calloc(R->num_nodes, sizeof(Router));
    assert(R->rl != NULL);
    for (int i = 0; i < R->num_nodes; ++i)
    {
        R->rl[i].id = i;
        R->rl[i].local_o.dst_id = -1;
        R->rl[i].channels = (Channel **)calloc(R->degree, sizeof(Channel *));
        assert(R->rl[i].channels != NULL);
        for (int j = 0; j < R->degree; ++j)
        {
            R->rl[i].channels[j] = (Channel *)calloc(ROUTE_SIM_ANYTOPO_VC, sizeof(Channel));
            assert(R->rl[i].channels[j] != NULL);
        }
        for (int j = 0; j < ROUTE_SIM_ANYTOPO_VC; ++j) // VC
        {
            for (int k = 0; k < R->degree; ++k)
            {
                R->rl[i].channels[k][j].buffer_o.dst_id = -1;
                if (G->adj_list[i][k] != INF)
                    R->rl[i].channels[k][j].credit = ROUTE_SIM_ANYTOPO_DEPTH;
            }
        }
    }
    for (int i = 0; i < G->num_nodes; ++i)
    {
        R->rl[i].neighbor = (Router **)calloc(G->degree, sizeof(Router *));
        assert(R->rl[i].neighbor != NULL);
        for (int j = 0; j < G->degree; ++j)
        {
            R->rl[i].neighbor[j] = &(R->rl[G->adj_list[i][j]]);
        }
    }
    int port;
    for (int i = 0; i < G->num_nodes; ++i)
    {
        for (int j = 0; j < G->degree; ++j)
        {
            if (G->adj_list[i][j] != INF && R->rl[i].channels[j][0].connect == NULL)
            {
                port = -1;
                for (int k = 0; k < G->degree; ++k)
                    if (G->adj_list[G->adj_list[i][j]][k] == i)
                    {
                        port = k;
                        break;
                    }
                assert(port != -1);
                for (int k = 0; k < ROUTE_SIM_ANYTOPO_VC; ++k)
                {
                    R->rl[i].channels[j][k].connect = &(R->rl[G->adj_list[i][j]].channels[port][k]);
                    R->rl[G->adj_list[i][j]].channels[port][k].connect = &(R->rl[i].channels[j][k]);
                }
            }
        }
    }
    // calloc non_local_congestion
    for (int i = 0; i < R->num_nodes; ++i) // node
    {
        for (int k = 0; k < ROUTE_SIM_ANYTOPO_VC; ++k)
        {
            R->rl[i].non_local_congestion[k] = (int *)calloc(R->degree, sizeof(int));
            assert(R->rl[i].non_local_congestion[k] != NULL);
            for (int j = 0; j < R->degree; ++j)
                if (G->adj_list[i][j] != INF)
                    R->rl[i].non_local_congestion[k][j] = ROUTE_SIM_ANYTOPO_DEPTH;
        }
    }
    return R;
}

void delete_router_list_traffic(RouterList *R)
{
    for (int i = 0; i < R->num_nodes; ++i)
    {
        if (R->rl[i].local_channel != NULL)
        {
            free(R->rl[i].local_channel);
            R->rl[i].local_channel = NULL;
        }
        R->rl[i].local_size = 0;
        R->rl[i].local_start = 0;
        R->rl[i].local_o.dst_id = -1;
        R->rl[i].recv_num = 0;
        R->rl[i].total_time = 0;
        for (int k = 0; k < R->degree; ++k)
        {
            for (int j = 0; j < ROUTE_SIM_ANYTOPO_VC; ++j)
            {
                R->rl[i].channels[k][j].buffer_o.dst_id = -1;
                R->rl[i].channels[k][j].size = 0;
                R->rl[i].channels[k][j].start = 0;
                R->rl[i].channels[k][j].end = 0;
                R->rl[i].channels[k][j].credit = ROUTE_SIM_ANYTOPO_DEPTH;
            }
        }
    }
}

void delete_router_list_lut(RouterList *R)
{
    for (int idx = 0; idx < ROUTE_SIM_ANYTOPO_VC; ++idx)
    {
        for (int i = 0; i < R->num_nodes; ++i)
        {
            if (R->rl[i].lut[idx] != NULL)
            {
                for (int j = 0; j < R->degree + 1; ++j)
                {
                    free(R->rl[i].lut[idx][j]);
                    R->rl[i].lut[idx][j] = NULL;
                }
                free(R->rl[i].lut[idx]);
                R->rl[i].lut[idx] = NULL;
            }
        }
    }
}

void delete_router_list(RouterList *R)
{
    for (int i = 0; i < R->num_nodes; ++i)
    {
        if (R->rl[i].local_channel != NULL)
        {
            free(R->rl[i].local_channel);
            R->rl[i].local_channel = NULL;
        }
        if (R->rl[i].channels != NULL)
        {
            for (int j = 0; j < R->degree; ++j)
            {
                if (R->rl[i].channels[j] != NULL)
                {

                    free(R->rl[i].channels[j]);
                    R->rl[i].channels[j] = NULL;
                }
            }
            free(R->rl[i].channels);
            R->rl[i].channels = NULL;
        }
    }
    delete_router_list_lut(R);
    free(R->rl);
    R->rl = NULL;
    free(R);
    R = NULL;
}

void init_router_list_traffic(RouterList *R, int **traffic_table)
{
    assert(traffic_table != NULL);
    int send_num;
    for (int i = 0; i < R->num_nodes; ++i)
    {
        send_num = 0;
        for (int j = 0; j < R->num_nodes; ++j)
        {
            if (traffic_table[i][j] > 0)
                send_num++;
        }
        R->rl[i].local_size = send_num;
        // printf("//// send from : %d to \t", i);
        if (send_num > 0)
        {
            R->rl[i].local_channel = (Packets *)calloc(send_num, sizeof(Packets));
            assert(R->rl[i].local_channel != NULL);
            for (int j = 0, k = 0; j < R->num_nodes; ++j)
            {
                if (traffic_table[i][j] > 0)
                {
                    // printf("(%d, %d) ", j, traffic_table[i][j]);
                    R->rl[i].local_channel[k].dst_id = j;
                    R->rl[i].local_channel[k].num = traffic_table[i][j];
                    ++k;
                }
            }
        }
        else
            R->rl[i].local_channel = NULL;
        // printf("\n");
    }
}

int **traffic_uniform_random(const RouterList *R, int p, int packet_num, unsigned int seed)
{
    srand(seed);
    int **traffic_table = (int **)calloc(R->num_nodes, sizeof(int *));
    assert(traffic_table != NULL);
    for (int i = 0; i < R->num_nodes; ++i)
    {
        traffic_table[i] = (int *)calloc(R->num_nodes, sizeof(int));
        assert(traffic_table[i] != NULL);
    }
    for (int i = 0; i < R->num_nodes; ++i)
    {
        for (int j = 0; j < R->num_nodes; ++j)
        {
            if ((rand() % 10000) < p)
                traffic_table[i][j] = packet_num;
        }
    }
    return traffic_table;
}

int **traffic_from_file(char *file_name, int num)
{
    // alloc
    int **traffic_table = (int **)calloc(num, sizeof(int *));
    assert(traffic_table != NULL);
    for (int i = 0; i < num; ++i)
    {
        traffic_table[i] = (int *)calloc(num, sizeof(int));
        assert(traffic_table[i] != NULL);
    }
    // read
    FILE *f = fopen(file_name, "r");
    assert(f != NULL);
    char *one_line = (char *)calloc(num * 15, sizeof(char));
    assert(one_line != NULL);
    char *token;
    int temp;
    for (int i = 0; i < num; ++i)
    {
        assert(fgets(one_line, sizeof(char) * num * 15, f) != NULL);
        token = strtok(one_line, ",");
        int j = 0;
        while (token != NULL)
        {
            temp = atoi(token);
            assert(temp >= 0);
            traffic_table[i][j] = temp;
            token = strtok(NULL, ",");
            ++j;
        }
    }
    fclose(f);
    f = NULL;
    free(one_line);
    one_line = NULL;
    token = NULL;
    return traffic_table;
}

void init_router_list_lut_shortest_path(RouterList *R, const GraphAdjList *G, int root_select, int (*turn_check)(int, int, int, int, int, int), unsigned char ***port_map)
{
    // Generate root & TC
    TreeCoordinate **tc = (TreeCoordinate **)calloc(ROUTE_SIM_ANYTOPO_VC, sizeof(TreeCoordinate *));
    assert(tc != NULL);
    int *root = (int *)calloc(ROUTE_SIM_ANYTOPO_VC, sizeof(int));
    assert(root != NULL);
    if (root_select == 0) // radnom select one root, and appply it to every VC
    {
        int best_root = rand() % G->num_nodes;
        for (int i = 0; i < ROUTE_SIM_ANYTOPO_VC; ++i)
        {
            root[i] = best_root;
            tc[i] = gen_tree_coordinate(G, root[i], 1);
        }
    }
    else if (root_select == 1) // random select one root for every VC independently
    {
        for (int i = 0; i < ROUTE_SIM_ANYTOPO_VC; ++i)
        {
            root[i] = rand() % G->num_nodes;
            tc[i] = gen_tree_coordinate(G, root[i], 1);
        }
    }
    else if (root_select == 2) // optimally select roots for every VC
    {
        root = root_selection(G, turn_check, ROUTE_SIM_ANYTOPO_VC);
        for (int i = 0; i < ROUTE_SIM_ANYTOPO_VC; ++i)
        {
            tc[i] = gen_tree_coordinate(G, root[i], 1);
            printf("\t best root is : %d", root[i]);
        }
        printf("\n");
        // tc = gen_tree_coordinate(G, best_root, 1);
    }
    else
        assert(0);

    // statistic path
    // int *path_sta = statistic_path(G, tc[0]);
    // printf("/*/*/*/*/*/ %d %d %d %d %d %d\n", path_sta[0], path_sta[1], path_sta[2], path_sta[3], path_sta[4], path_sta[5]);

    // generate LUT
    delete_router_list_lut(R);
    // find shortest path
    for (int i = 0; i < R->num_nodes; ++i)
    {
        memset((*port_map)[i], 255, R->num_nodes * sizeof(unsigned char));
        for (unsigned char j = 0; j < G->degree; ++j)
        {
            if (G->adj_list[i][j] != INF)
            {
                (*port_map)[i][G->adj_list[i][j]] = j; // dst[i][p] = j; means: node i's j-th port is connected to node p;
            }
        }
    }
    trun_model_bfs_new(G, tc, R, turn_check);
    for (int k = 0; k < ROUTE_SIM_ANYTOPO_VC; ++k)
    {
        // check lut
        int false_path_num = 0;
        for (int core = 0; core < R->num_nodes; ++core)
        {
            for (int i = 0; i < R->degree + 1; ++i)
            {
                for (int j = 0; j < R->num_nodes; ++j)
                {
                    if (R->rl[core].lut[k][i][j] == 255)
                    {
                        assert(R->rl[core].self_dst[k][i][j] == INF);
                        false_path_num++;
                    }
                    else
                    {
                        assert(R->rl[core].self_dst[k][i][j] != INF);
                    }
                    if (i == R->degree)
                    {
                        assert(R->rl[core].lut[k][i][j] != 255);
                        assert(R->rl[core].self_dst[k][i][j] != INF);
                    }
                }
            }
        }
        double false_dist_rate = (double)false_path_num / (double)(G->degree * G->num_nodes * G->num_nodes);
        printf("False dist rate is %lf\n", false_dist_rate);
        // check lut & self_dst by path
        int cur_dst = 0;
        int cur_node = 0;
        int out_port = 0;
        int next_node = 0;
        int in_port = 0;
        for (int i = 0; i < R->num_nodes; ++i) // src node
        {
            for (int j = 0; j < R->num_nodes; ++j) // dst node
            {
                if (i == j)
                {
                    assert(R->rl[i].lut[k][R->degree][j] == 254);
                    assert(R->rl[i].self_dst[k][R->degree][j] == 0);
                }
                else
                {
                    cur_node = i;
                    cur_dst = R->rl[cur_node].self_dst[k][R->degree][j];
                    out_port = R->rl[cur_node].lut[k][R->degree][j];
                    next_node = G->adj_list[cur_node][out_port];
                    in_port = (*port_map)[next_node][cur_node];
                    while (next_node != j)
                    {
                        assert(R->rl[next_node].self_dst[k][in_port][j] == cur_dst - 1);
                        cur_dst--;
                        cur_node = next_node;
                        out_port = R->rl[cur_node].lut[k][in_port][j];
                        next_node = G->adj_list[cur_node][out_port];
                        in_port = (*port_map)[next_node][cur_node];
                    }
                    assert(R->rl[j].lut[k][in_port][j] == 254);
                    assert(R->rl[j].self_dst[k][in_port][j] == 0);
                }
            }
        }
        // generate direction_viable
        unsigned int *neighbors = NULL;
        for (int i = 0; i < R->num_nodes; ++i) // node
        {
            R->rl[i].direction_viable[k] = (unsigned char *)calloc(R->degree, sizeof(unsigned char));
            assert(R->rl[i].direction_viable[k] != NULL);
            for (int j = 0; j < R->degree; ++j) // input port
            {
                for (int p = 0; p < R->degree; ++p) // output port
                {
                    neighbors = G->adj_list[i];
                    if (j == p || neighbors[j] == INF || neighbors[p] == INF)
                        continue;
                    if ((*turn_check)(tc[k]->x[neighbors[j]], tc[k]->y[neighbors[j]], tc[k]->x[i], tc[k]->y[i], tc[k]->x[neighbors[p]], tc[k]->y[neighbors[p]]) != 0)
                        R->rl[i].direction_viable[k][j] |= (1 << p); // bitmap
                }
            }
        }
        neighbors = NULL;
    }
    for (int i = 0; i < ROUTE_SIM_ANYTOPO_VC; ++i)
        delete_tc(tc[i]);
    free(tc);
    free(root);
    printf("--------------------- End of LUT INIT -------------------------\n");
}

double avg_packet_delivery_latency(const RouterList *R)
{
    unsigned long long total_time = 0;
    int total_packet_num = 0;
    for (int i = 0; i < R->num_nodes; ++i)
    {
        total_time += R->rl[i].total_time;
        assert(total_time < (1ULL << 62)); // in case of overflow
        total_packet_num += R->rl[i].recv_num;
    }
    return (double)((long double)total_time / (long double)total_packet_num);
}

void router_arbit(RouterList *R, int *router_empty, int *router_busy, unsigned int sim_time, const int path_diversity_mode, const int load_balance_mode, unsigned char **port_map, unsigned int **path_log)
{
    int p_dst_id, p_src_id;
    unsigned int p_start_time;
    int p_ttl;
    Router *r = NULL;
    unsigned char dst_port;
    int dst_vc;
    // load balance
    int viable_output_num;
    int *viable_output_list = (int *)calloc(R->degree, sizeof(int));
    assert(viable_output_list != NULL);
    for (int node = 0; node < R->num_nodes; ++node) // every node
    {
        r = &(R->rl[node]);
        for (int i = 0; i < R->degree; ++i) // every input direction's channel
        {
            for (int j = 0; j < ROUTE_SIM_ANYTOPO_VC; ++j) // every VC
            {
                if (&(r->channels[i][j]) != NULL && r->channels[i][j].size > 0) // channel's buffer_i not empty
                {
                    *router_empty = 0; // not all packets have arrived destination
                    p_src_id = r->channels[i][j].buffer_i[r->channels[i][j].start].src_id;
                    p_dst_id = r->channels[i][j].buffer_i[r->channels[i][j].start].dst_id;
                    p_start_time = r->channels[i][j].buffer_i[r->channels[i][j].start].start_time;
                    p_ttl = r->channels[i][j].buffer_i[r->channels[i][j].start].ttl;
                    // Path Diversity  &  Load Balance //
                    // determine the output port by "path_diversity_mode" & "load_balance_mode"
                    if (p_dst_id != node)
                    {
                        if (path_diversity_mode == -1) // **** no diversity ****
                            dst_port = r->lut[j][i][p_dst_id];
                        else if (path_diversity_mode == 0) //  ****minimal diversity: only minimal-path can be choosed ****
                        {
                            int shortest_dst = r->self_dst[j][i][p_dst_id];
                            if (load_balance_mode == 0) // equal load
                            {
                                viable_output_num = 0;
                                for (int output_idx = 0; output_idx < R->degree; ++output_idx)
                                {
                                    if (r->neighbor[output_idx]->self_dst[j][port_map[r->neighbor[output_idx]->id][r->id]][p_dst_id] == shortest_dst - 1)
                                        viable_output_list[viable_output_num++] = output_idx;
                                }
                                dst_port = viable_output_list[rand() % viable_output_num];
                            }
                            else if (load_balance_mode == 1) // local congestion aware
                            {
                                int biggest_credit = 0;
                                viable_output_num = 0;
                                // find the output channels, whoes credit is the biggest
                                for (int output_idx = 0; output_idx < R->degree; ++output_idx)
                                {
                                    if (r->neighbor[output_idx]->self_dst[j][port_map[r->neighbor[output_idx]->id][r->id]][p_dst_id] == shortest_dst - 1)
                                    {
                                        if (biggest_credit < r->channels[output_idx][j].credit)
                                        {
                                            viable_output_num = 0;
                                            viable_output_list[viable_output_num++] = output_idx;
                                            biggest_credit = r->channels[output_idx][j].credit;
                                        }
                                        else if (biggest_credit == r->channels[output_idx][j].credit && biggest_credit != 0)
                                        {
                                            viable_output_list[viable_output_num++] = output_idx;
                                        }
                                    }
                                }
                                if (viable_output_num == 0)
                                    dst_port = r->lut[j][i][p_dst_id];
                                else
                                    dst_port = viable_output_list[rand() % viable_output_num];
                            }
                            else if (load_balance_mode == 2) // non-local congestion aware
                            {
                                int biggest_credit = 0;
                                viable_output_num = 0;
                                // find the output channels, whoes credit is the biggest
                                for (int output_idx = 0; output_idx < R->degree; ++output_idx)
                                {
                                    if (r->neighbor[output_idx]->self_dst[j][port_map[r->neighbor[output_idx]->id][r->id]][p_dst_id] == shortest_dst - 1 && r->channels[output_idx][j].buffer_o.dst_id == -1) // only select free channel
                                    {
                                        // assert(r->non_local_congestion[j][output_idx] == r->channels[output_idx][j].credit);
                                        if (biggest_credit < r->non_local_congestion[j][output_idx])
                                        {
                                            viable_output_num = 0;
                                            viable_output_list[viable_output_num++] = output_idx;
                                            biggest_credit = r->non_local_congestion[j][output_idx];
                                        }
                                        else if (biggest_credit == r->non_local_congestion[j][output_idx] && biggest_credit != 0)
                                        {
                                            viable_output_list[viable_output_num++] = output_idx;
                                        }
                                    }
                                }
                                if (viable_output_num == 0)
                                    dst_port = r->lut[j][i][p_dst_id];
                                else
                                    dst_port = viable_output_list[rand() % viable_output_num];
                            }
                            else
                                assert(0);
                        }
                        else // **** non-minimal diversity ****
                        {
                            int shortest_dst = r->self_dst[j][i][p_dst_id];
                            if (load_balance_mode == 0) // equal load
                            {
                                assert(0);
                            }
                            else if (load_balance_mode == 1) // local congestion aware
                            {
                                int non_minimal_biggest_credit = 0;
                                int non_minimal_best_port = -1;
                                int minimal_biggest_credit = 0; // minimal path
                                viable_output_num = 0;
                                // find the minimal path output channels, whoes credit is the biggest
                                for (int output_idx = 0; output_idx < R->degree; ++output_idx)
                                {
                                    if (r->neighbor[output_idx]->self_dst[j][port_map[r->neighbor[output_idx]->id][r->id]][p_dst_id] == shortest_dst - 1)
                                    {
                                        if (minimal_biggest_credit < r->channels[output_idx][j].credit)
                                        {
                                            viable_output_num = 0;
                                            viable_output_list[viable_output_num++] = output_idx;
                                            minimal_biggest_credit = r->channels[output_idx][j].credit;
                                        }
                                        else if (minimal_biggest_credit == r->channels[output_idx][j].credit && minimal_biggest_credit != 0)
                                        {
                                            viable_output_list[viable_output_num++] = output_idx;
                                        }
                                    }
                                    else
                                    {
                                        if (r->neighbor[output_idx]->self_dst[j][port_map[r->neighbor[output_idx]->id][r->id]][p_dst_id] != INF)
                                        {
                                            if (non_minimal_biggest_credit < r->channels[output_idx][j].credit)
                                            {
                                                non_minimal_biggest_credit = r->channels[output_idx][j].credit;
                                                non_minimal_best_port = output_idx;
                                            }
                                            else if (non_minimal_biggest_credit == r->channels[output_idx][j].credit && non_minimal_biggest_credit != 0)
                                            {
                                                if (rand() % 2 == 0)
                                                {
                                                    non_minimal_biggest_credit = r->channels[output_idx][j].credit;
                                                    non_minimal_best_port = output_idx;
                                                }
                                            }
                                        }
                                    }
                                }
                                if (p_ttl > 0)
                                {
                                    // if (non_minimal_biggest_credit > minimal_biggest_credit) // choose non-minimal path
                                    if (non_minimal_biggest_credit == ROUTE_SIM_ANYTOPO_DEPTH && minimal_biggest_credit == 0)
                                    {
                                        dst_port = non_minimal_best_port;
                                        p_ttl--;
                                    }
                                    else
                                    { // choose minimal path
                                        if (viable_output_num == 0)
                                            dst_port = r->lut[j][i][p_dst_id];
                                        else
                                            dst_port = viable_output_list[rand() % viable_output_num];
                                    }
                                }
                                else
                                { // choose minimal path
                                    if (viable_output_num == 0)
                                        dst_port = r->lut[j][i][p_dst_id];
                                    else
                                        dst_port = viable_output_list[rand() % viable_output_num];
                                }
                            }
                            else if (load_balance_mode == 2) // non-local congestion aware
                            {
                                int non_minimal_biggest_credit = 0;
                                int non_minimal_best_port = -1;
                                int minimal_biggest_credit = 0; // minimal path
                                viable_output_num = 0;
                                // find the minimal path output channels, whoes credit is the biggest
                                for (int output_idx = 0; output_idx < R->degree; ++output_idx)
                                {
                                    if (r->neighbor[output_idx]->self_dst[j][port_map[r->neighbor[output_idx]->id][r->id]][p_dst_id] == shortest_dst - 1 && r->channels[output_idx][j].buffer_o.dst_id == -1) // only select free channel
                                    {
                                        if (minimal_biggest_credit < r->channels[output_idx][j].credit)
                                        {
                                            viable_output_num = 0;
                                            viable_output_list[viable_output_num++] = output_idx;
                                            minimal_biggest_credit = r->channels[output_idx][j].credit;
                                        }
                                        else if (minimal_biggest_credit == r->channels[output_idx][j].credit && minimal_biggest_credit != 0)
                                        {
                                            viable_output_list[viable_output_num++] = output_idx;
                                        }
                                    }
                                    else
                                    {
                                        if (r->neighbor[output_idx]->self_dst[j][port_map[r->neighbor[output_idx]->id][r->id]][p_dst_id] != INF && r->channels[output_idx][j].buffer_o.dst_id == -1) // only select free channel
                                        {
                                            if (non_minimal_biggest_credit < r->non_local_congestion[j][output_idx])
                                            {
                                                non_minimal_biggest_credit = r->non_local_congestion[j][output_idx];
                                                non_minimal_best_port = output_idx;
                                            }
                                            else if (non_minimal_biggest_credit == r->non_local_congestion[j][output_idx] && non_minimal_biggest_credit != 0)
                                            {
                                                if (rand() % 2 == 0)
                                                {
                                                    non_minimal_biggest_credit = r->non_local_congestion[j][output_idx];
                                                    non_minimal_best_port = output_idx;
                                                }
                                            }
                                        }
                                    }
                                }
                                if (p_ttl > 0)
                                {
                                    // if (non_minimal_biggest_credit >= (R->degree * R->degree * R->degree * ROUTE_SIM_ANYTOPO_DEPTH) && minimal_biggest_credit == 0)
                                    if (non_minimal_best_port != -1 && r->channels[non_minimal_best_port][j].credit == ROUTE_SIM_ANYTOPO_DEPTH && minimal_biggest_credit == 0)
                                    {
                                        dst_port = non_minimal_best_port;
                                        p_ttl--;
                                    }
                                    else
                                    { // choose minimal path
                                        if (viable_output_num == 0)
                                            dst_port = r->lut[j][i][p_dst_id];
                                        else
                                            dst_port = viable_output_list[rand() % viable_output_num];
                                    }
                                }
                                else
                                { // choose minimal path
                                    if (viable_output_num == 0)
                                        dst_port = r->lut[j][i][p_dst_id];
                                    else
                                        dst_port = viable_output_list[rand() % viable_output_num];
                                }
                            }
                            else
                                assert(0);
                        }
                    }
                    // printf("---- packet in NODE %d, channel %d, DST is %d, dst_port %d \t", node, i, p_dst_id, dst_port);
                    if (p_dst_id == node && r->local_o.dst_id == -1) // move to local
                    {
                        // printf("\t ** moved !\n");
                        r->local_o.src_id = p_src_id;
                        r->local_o.dst_id = p_dst_id;
                        r->local_o.start_time = p_start_time;
                        r->local_o.ttl = p_ttl;
                        r->channels[i][j].size -= 1;
                        r->channels[i][j].start = (r->channels[i][j].start + 1) % ROUTE_SIM_ANYTOPO_DEPTH; // pop this packet
                        *router_busy = 1;                                                                  // at least onr packet is moving
                        r->recv_num++;
                        // r->total_time += sim_time;
                        r->total_time += sim_time - p_start_time + 1;
                        assert(r->total_time < (1ULL << 62)); // in case of overflow
                    }
                    else if (p_dst_id != node && r->channels[dst_port][j].buffer_o.dst_id == -1 && r->channels[dst_port][j].credit > 0) // move to dst_port
                    {
                        // printf("\t ** moved !\n");
                        r->channels[dst_port][j].buffer_o.src_id = p_src_id; // assign to dst_port's buffer_o
                        r->channels[dst_port][j].buffer_o.dst_id = p_dst_id; // assign to dst_port's buffer_o
                        r->channels[dst_port][j].buffer_o.start_time = p_start_time;
                        r->channels[dst_port][j].buffer_o.ttl = p_ttl;
                        r->channels[dst_port][j].credit -= 1;
                        r->channels[i][j].size -= 1;
                        r->channels[i][j].start = (r->channels[i][j].start + 1) % ROUTE_SIM_ANYTOPO_DEPTH; // pop this packet
                        *router_busy = 1;                                                                  // at least one packet is moving
                        path_log[r->id][dst_port] += 1;                                                    // record the number of delivered packets for every path
                    }
                    // else
                    //     printf("\n");
                }
            }
        }
        if (r->local_start < r->local_size) // local channel
        {
            *router_empty = 0;
            p_src_id = node;
            p_dst_id = r->local_channel[r->local_start].dst_id;
            dst_vc = rand() % ROUTE_SIM_ANYTOPO_VC; // random send to each VC
            dst_port = r->lut[dst_vc][R->degree][p_dst_id];
            if (p_dst_id == node && r->local_o.dst_id == -1) // move to local
            {
                r->local_o.src_id = p_src_id;
                r->local_o.dst_id = p_dst_id;
                // r->local_o.start_time = sim_time;
                // r->local_o.ttl = 20;  // TTL
                r->local_channel[r->local_start].num--; // pop this packet
                if (r->local_channel[r->local_start].num == 0)
                    r->local_start++;
                *router_busy = 1;
                r->recv_num++;
                // r->total_time += sim_time;
                r->total_time += 1;
                assert(r->total_time < (1ULL << 62)); // in case of overflow
            }
            else if (p_dst_id != node && r->channels[dst_port][dst_vc].buffer_o.dst_id == -1 && r->channels[dst_port][dst_vc].credit > 0) // move to dst_port
            {
                r->channels[dst_port][dst_vc].buffer_o.src_id = p_src_id; // assign to dst_port's buffer_o
                r->channels[dst_port][dst_vc].buffer_o.dst_id = p_dst_id; // assign to dst_port's buffer_o
                r->channels[dst_port][dst_vc].buffer_o.start_time = sim_time;
                r->channels[dst_port][dst_vc].buffer_o.ttl = ROUTE_SIM_TTL; // TTL
                r->channels[dst_port][dst_vc].credit -= 1;
                r->local_channel[r->local_start].num--; // pop this packet
                if (r->local_channel[r->local_start].num == 0)
                    r->local_start++;
                *router_busy = 1;
                path_log[r->id][dst_port] += 1; // record the number of delivered packets for every path
            }
        }
    }
    free(viable_output_list);
}

void router_trans(RouterList *R, const int path_diversity_mode, const int load_balance_mode, unsigned char **port_map)
{
    Router *r = NULL;
    for (int node = 0; node < R->num_nodes; ++node) // every node
    {
        r = &(R->rl[node]);
        for (int i = 0; i < R->degree; ++i)
        {
            for (int j = 0; j < ROUTE_SIM_ANYTOPO_VC; ++j)
            {
                if (&(r->channels[i][j]) != NULL && r->channels[i][j].connect != NULL)
                {
                    if (r->channels[i][j].connect->buffer_o.dst_id != -1) // move packet from buffer_o to buffer_i
                    {
                        assert(r->channels[i][j].size < ROUTE_SIM_ANYTOPO_DEPTH);
                        r->channels[i][j].buffer_i[r->channels[i][j].end].src_id = r->channels[i][j].connect->buffer_o.src_id;
                        r->channels[i][j].buffer_i[r->channels[i][j].end].dst_id = r->channels[i][j].connect->buffer_o.dst_id;
                        r->channels[i][j].buffer_i[r->channels[i][j].end].start_time = r->channels[i][j].connect->buffer_o.start_time;
                        r->channels[i][j].buffer_i[r->channels[i][j].end].ttl = r->channels[i][j].connect->buffer_o.ttl;
                        r->channels[i][j].end = (r->channels[i][j].end + 1) % ROUTE_SIM_ANYTOPO_DEPTH;
                        r->channels[i][j].size++;
                        r->channels[i][j].connect->buffer_o.dst_id = -1; // clear
                    }
                    r->channels[i][j].connect->credit = ROUTE_SIM_ANYTOPO_DEPTH - r->channels[i][j].size; // update credit
                }
            }
        }
        r->local_o.dst_id = -1; // clear local
    }
    if (path_diversity_mode != -1 && load_balance_mode == 2)
    {
        // update non-local congestion infomation
        Router *l1_node, *l2_node, *l3_node;
        for (int i = 0; i < R->degree; ++i) // output port
        {
            l1_node = r->neighbor[i]; // next node
            for (int k = 0; k < ROUTE_SIM_ANYTOPO_VC; ++k)
            {
                int l1_rate = 0;
                r->non_local_congestion[k][i] = 0;
                for (int j = 0; j < R->degree; ++j)
                {
                    if (((l1_node->direction_viable[k][port_map[l1_node->id][r->id]] >> j) & 0x1) == 1)
                    {
                        l2_node = l1_node->neighbor[j];
                        if (l2_node == r || l2_node == NULL)
                            continue;
                        int l2_rate = 0;
                        int l2_credit_temp = 0;
                        for (int m = 0; m < R->degree; ++m)
                        {
                            if (((l2_node->direction_viable[k][port_map[l2_node->id][l1_node->id]] >> m) & 0x1) == 1)
                            {
                                l3_node = l2_node->neighbor[m];
                                if (l3_node == r || l3_node == l1_node || l3_node == NULL)
                                    continue;
                                l2_credit_temp += l2_node->channels[m][k].credit;
                                l2_rate++;
                            }
                        }
                        if (l1_node->channels[j][k].credit != 0)
                        {
                            r->non_local_congestion[k][i] += l2_credit_temp;
                            r->non_local_congestion[k][i] += l1_node->channels[j][k].credit * l2_rate;
                            l1_rate += l2_rate;
                        }
                    }
                }
                if (r->channels[i][k].credit == 0)
                {
                    r->non_local_congestion[k][i] = 0;
                }
                else
                {
                    r->non_local_congestion[k][i] += r->channels[i][k].credit * l1_rate;
                }
            }
        }
    }
}

void path_log_initial(const GraphAdjList *G, unsigned int ***path_log)
{
    (*path_log) = (unsigned int **)calloc(G->num_nodes, sizeof(unsigned int *));
    assert((*path_log) != NULL);
    for (int i = 0; i < G->num_nodes; ++i)
    {
        (*path_log)[i] = (unsigned int *)calloc(G->degree, sizeof(unsigned int));
        assert((*path_log)[i] != NULL);
    }
}

void path_log_delete(const GraphAdjList *G, unsigned int **path_log)
{
    for (int i = 0; i < G->num_nodes; ++i)
    {
        free(path_log[i]);
        path_log[i] = NULL;
    }
    free(path_log);
    path_log = NULL;
}

void route_sim_anytopo(char *file_name, int traffic_mode, int route_lut_mode, unsigned int seed, int ppp, int packets_num, int root_select, int path_diversity_mode, int load_balance_mode, char *path_log_file_name, char *traffic_name, int traffic_num)
{
    printf("\n--------------------- Start of Rounte_Sim_Anytopo -------------------------\n");
    struct timeval time_start, time_end;

    // print parameters
    printf("**** Graph name is : %s\n", file_name);
    printf("**** Traffic Mode is : %d\n", traffic_mode);                                 // -1:read from file; 0-uniform random
    printf("     seed = %u; ppp = %d/10000; packet_num = %d\n", seed, ppp, packets_num); // 0-uniform random
    printf("     traffic LUT file name is: %s\n", traffic_name);
    printf("**** Route Lut Mode is : %d\n", route_lut_mode);           // 0-L-turn; 1-Tree-turn; 2-new-turn
    printf("**** Root selection mode is : %d\n", root_select);         // 0-fix; 1-random; 2-optimal select
    printf("**** Path Diversity mode is : %d\n", path_diversity_mode); // -1-no diversity; 0-minimal-diversity; 1-non-minimal diversity
    printf("**** Load Balance mode is : %d\n", load_balance_mode);     // 0-equal load; 1-local congestion aware; 2-non-local congestion aware
    printf("**** Num of VC is : %d\n", ROUTE_SIM_ANYTOPO_VC);
    printf("**** Num of LUT is : %d\n", ROUTE_SIM_ANYTOPO_VC);
    printf("**** Depth of Channel is : %d\n", ROUTE_SIM_ANYTOPO_DEPTH);
    printf("**** Path-log file name is : %s\n", path_log_file_name);

    GraphAdjList *G = construct_graph_adjlist_from_file(file_name);
    // GraphAdjList *G = jellyfish_gen(1100, 6, 0);
    // GraphAdjList *G = torus_gen(20, 2);
    // print_graph_adjlist(G);
    check_graph_adjlist(G);
    RouterList *R = construct_router_list_from_graphadjlist(G);

    //
    unsigned char **port_map = (unsigned char **)calloc(G->num_nodes, sizeof(unsigned char *));
    assert(port_map != NULL);
    for (int i = 0; i < R->num_nodes; ++i)
    {
        port_map[i] = (unsigned char *)calloc(R->num_nodes, sizeof(unsigned char));
        assert(port_map[i] != NULL);
    }

    // Traffic Mode
    int **traffic_table = NULL;
    if (traffic_mode == -1)
        traffic_table = traffic_from_file(traffic_name, traffic_num);
    else if (traffic_mode == 0)
        traffic_table = traffic_uniform_random(R, ppp, packets_num, 0);
    else
        assert(0);

    init_router_list_traffic(R, traffic_table);

    // Route Lut Mode
    if (route_lut_mode == 0)
        init_router_list_lut_shortest_path(R, G, root_select, check_l_trun, &port_map);
    else if (route_lut_mode == 1)
        init_router_list_lut_shortest_path(R, G, root_select, check_tree_trun, &port_map);
    else if (route_lut_mode == 2)
        init_router_list_lut_shortest_path(R, G, root_select, check_octo_trun, &port_map);
    else
        assert(0);

    // path log
    unsigned int **path_log = NULL;
    path_log_initial(G, &path_log);
    // start sim
    unsigned sim_time = 0;
    int dead_lock_cnt = 0;
    int router_empty = 1;
    int router_busy = 0;
    gettimeofday(&time_start, NULL);
    printf("Start ...\n");
    while (1)
    {
        router_arbit(R, &router_empty, &router_busy, sim_time, path_diversity_mode, load_balance_mode, port_map, path_log);
        router_trans(R, path_diversity_mode, load_balance_mode, port_map);
        // if (sim_time % 10000 == 0)
        //     printf("-- sim_time : %d\n", sim_time);
        if (router_busy == 0)
        {
            if (router_empty == 1)
            {
                printf("Total cycle : %u\n", sim_time);
                printf("End ...\n");
                break;
            }
            else if (dead_lock_cnt > 10)
            {
                printf("Total cycle : %u\n", sim_time);
                printf("****** ERROR ! Dead Lock ...\n");
                break;
            }
            else
                dead_lock_cnt += 1;
        }
        else
            dead_lock_cnt = 0;
        sim_time += 1;
        router_empty = 1;
        router_busy = 0;
    }
    gettimeofday(&time_end, NULL);
    double total_time = ((double)(time_end.tv_sec) + (double)(time_end.tv_usec) / 1000000.0) - ((double)(time_start.tv_sec) + (double)(time_start.tv_usec) / 1000000.0);
    printf("Using time = %f s\n", total_time);
    printf("**** Sim Done !\n");

    printf("Avg Packet Delivery Latency is : %lf\n", avg_packet_delivery_latency(R));

    write_path_log_to_file(G, path_log, path_log_file_name);

    // FREE
    for (int i = 0; i < R->num_nodes; ++i)
        free(port_map[i]);
    free(port_map);
    path_log_delete(G, path_log);
}