#ifndef __ROUTE_SIM_ANYTOPO_H__
#define __ROUTE_SIM_ANYTOPO_H__ 1

#include "topology.h"

#define ROUTE_SIM_ANYTOPO_DEPTH 3
#define ROUTE_SIM_ANYTOPO_VC 2
#define ROUTE_SIM_TTL 20

typedef struct
{
    int dst_id;
    int num;
} Packets;

typedef struct
{
    int src_id;
    int dst_id;
    unsigned int start_time;
    int ttl;
} Packet; // -1 for invalid packet

typedef struct channel Channel;
struct channel
{
    Packet buffer_i[ROUTE_SIM_ANYTOPO_DEPTH];
    Packet buffer_o;
    Channel *connect; // connect to another Router's Channel
    unsigned char size;
    unsigned char start;
    unsigned char end;
    unsigned char credit;
};

typedef struct router Router;
struct router
{
    int id;
    // local channel
    Packets *local_channel;
    int local_size;
    int local_start;
    Packet local_o;
    // other channels
    Channel **channels; // [degree][VC]
    // time statistic
    int recv_num;
    unsigned long long total_time;
    // LUT
    unsigned char **lut[ROUTE_SIM_ANYTOPO_VC]; // [VC][input_degree + 1][num_nodes(dst)]
    // for each input channel, the viable output direction of output
    unsigned char *direction_viable[ROUTE_SIM_ANYTOPO_VC];
    // the shortest distance to every node from this node
    unsigned int **self_dst[ROUTE_SIM_ANYTOPO_VC];
    Router **neighbor;
    // non-local congestion infomation
    int *non_local_congestion[ROUTE_SIM_ANYTOPO_VC];
};

typedef struct
{
    int num_nodes;
    int degree;
    Router *rl;
} RouterList;

RouterList *construct_router_list_from_graphadjlist(GraphAdjList *G);

void delete_router_list_traffic(RouterList *R);

void delete_router_list_lut(RouterList *R);

void delete_router_list(RouterList *R);

void init_router_list_traffic(RouterList *R, int **traffic_table);

int **traffic_uniform_random(const RouterList *R, int p, int packet_num, unsigned int seed);

void route_sim_anytopo(char *file_name, int traffic_mode, int route_lut_mode, unsigned int seed, int ppp, int packets_num, int root_select, int path_diversity_mode, int load_balance_mode, char *path_log_file_name, char *traffic_name, int traffic_num);

double avg_packet_delivery_latency(const RouterList *R);

#endif