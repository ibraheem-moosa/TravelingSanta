#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>

#define N 197770 
#define GRID_LEN 128
#define MAXX 5100
#define MAXY 3400
#define MAX_DIST 7000
#define STARTING_POS 121572

struct pos
{
    int i;
    int j;
};

int da_x[] = { 0, -1, -1, 0, 1, 1, -1,  1};
int da_y[] = {-1, -1,  0, 1, 1, 0,  1, -1};

void get_chebyshev_neighbors(struct pos ij, int mx, int my, int* n, struct pos* ret)
{
    int i = ij.i;
    int j = ij.j;
    *n = 8;
    int ci = 0;
    for(int k = 0; k < 8; k++)
    {
        int di = i + da_x[k];
        int dj = j + da_y[k];
        if(di < 0 || di > mx) continue;
        if(dj < 0 || dj > my) continue;
        ret[ci].i = di;
        ret[ci].j = dj;
        ci++;
    }
    *n = ci;
}

struct pos get_pos_from_xy(double x, double y)
{
    struct pos ret;
    ret.i = (int)floor(x / MAXX);
    ret.j = (int)floor(y / MAXY);
    return ret;
}

static char is_prime[N];

int get_all_primes(int below, int** out)
{
    int prime_count = 0;
    for(int i = 0; i < below; i++) is_prime[i] = 1;
    is_prime[0] = 0;
    is_prime[1] = 0;
    for(int i = 2; i < below; i++)
    {
        if(is_prime[i])
        {
            prime_count++;
            for(int q = 2 * i; q < below; q += i)
                is_prime[q] = 0;
        }
    }
    int* ret = malloc(prime_count * sizeof(int));
    int current_prime_index = 0;
    for(int i = 0; i < below; i++)
    {
        if(is_prime[i])
            ret[current_prime_index++] = i;
    }
    *out = ret;
    return prime_count;
}

void random_permutattion(int* array, int n)
{
    for(int i = 0; i < n - 1; i++)
    {
        int j = i + rand() % (n - i);
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
} 

double euclidean_distance(double x1, double y1, double x2, double y2)
{
    double dx = x1 - x2;
    double dy = y1 - y2;
    return sqrt(dx * dx + dy * dy); 
}


struct city
{
    int id;
    double x;
    double y;
};

struct city cities[N];

void load_cities(const char* file_name)
{
    FILE* input_file = fopen(file_name, "r");
    char* buf = NULL;
    size_t buf_len = 0;
    while(1)
    {
        ssize_t read_count = 
            getline(&buf, &buf_len, input_file);
        if(read_count == -1)
            break;
        int city_id;
        double x, y;
        sscanf(buf, "%d,%lf,%lf", &city_id, &x, &y);
        cities[city_id].id = city_id;
        cities[city_id].x = x;
        cities[city_id].y = y;
    }
    free(buf);
    fclose(input_file);
}

double dist_between_id(int id1, int id2)
{
    double x1 = cities[id1].x;
    double y1 = cities[id1].y;
    double x2 = cities[id2].x;
    double y2 = cities[id2].y;
    return euclidean_distance(x1, y1, x2, y2);
}


struct bucket_node
{
    int id;
    struct bucket_node* next;
};

struct bucket_node* init_bucket_node(int id)
{
    struct bucket_node* ret = malloc(sizeof(struct bucket_node));
    ret->id = id;
    ret->next = NULL;
    return ret;
}

void prepend_bucket_node(struct bucket_node* bn, struct bucket_node* before)
{
    bn->next = before;
}

void delete_bucket_node(struct bucket_node* bn)
{
    if(bn) {
        delete_bucket_node(bn->next);
        free(bn);
    }
}

struct bucket_grid_2d
{
    int grid_len;
    struct bucket_node** buckets;
};

struct bucket_grid_2d init_bucket_grid_2d(int grid_len)
{
    struct bucket_grid_2d ret;
    ret.grid_len = grid_len;
    ret.buckets = calloc(grid_len * grid_len, sizeof(struct bucket_node*));
    return ret;
}

struct bucket_node* get_at_bucket_grid_2d(struct bucket_grid_2d bg, int i, int j)
{
    int pos = bg.grid_len * i + j;
    return bg.buckets[pos];
}

void insert_at_bucket_grid_2d(struct bucket_grid_2d bg, int i, int j, int id)
{
    int pos = bg.grid_len * i + j;
    struct bucket_node* bn = init_bucket_node(id);
    struct bucket_node* before = bg.buckets[pos];
    if(before) {
        prepend_bucket_node(bn, before);
    }
    bg.buckets[pos] = bn;
}

void remove_from_bucket_grid_2d(struct bucket_grid_2d bg, int i, int j, int id)
{
    int pos = bg.grid_len * i + j;
    struct bucket_node* bn = bg.buckets[pos];
    if(bn->id == id)
    {
        bg.buckets[pos] = bn->next;
        free(bn);
        return;
    }
    struct bucket_node* prev = bn;
    bn = bn->next;
    while(bn)
    {
        if(bn->id == id)
        {
            prev->next = bn->next;
            free(bn);
            break;
        }
        prev = bn;
        bn = bn->next;
    }
}

void delete_bucket_grid_2d(struct bucket_grid_2d bg)
{
    int n = bg.grid_len;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            struct bucket_node* bn = get_at_bucket_grid_2d(bg, i, j);
            delete_bucket_node(bn);
        }
    }
}

static struct bucket_grid_2d bucket_grid;

struct solution
{
    int n;
    int* ids;
};

struct solution load_solution(char* fname, int n)
{
    struct solution ret;
    ret.n = n;
    ret.ids = malloc(n * sizeof(int));
    FILE* f = fopen(fname, "r");
    char* buf = NULL;
    size_t buf_len = 0;
    int current_index = 0;
    //ignore first line
    getline(&buf, &buf_len, f);
    while(1)
    {
        ssize_t read_count = 
            getline(&buf, &buf_len, f);
        if(read_count == -1)
            break;
        sscanf(buf, "%d", &ret.ids[current_index++]);
    }
    free(buf);
    fclose(f);
    return ret;
}


int get_nearest_neighbor(int id, char* visited, struct bucket_grid_2d bg, char prime_preferred)
{
    struct pos ij = get_pos_from_xy(cities[id].x, cities[id].y);
    int num_neighbors;
    struct pos neighbors[8];
    get_chebyshev_neighbors(ij, MAXX, MAXY, &num_neighbors, neighbors);
    double min_dist = MAX_DIST;
    int min_id = 0;
    int ci = 0;
    do 
    {
        struct bucket_node* bn = get_at_bucket_grid_2d(bg, ij.i, ij.j);
        while(bn)
        {
            if(!visited[bn->id])
            {
                double dist = dist_between_id(id, bn->id);
                if(prime_preferred && !is_prime[bn->id])
                    dist *= 1.1;
                if (min_dist > dist)
                {
                    min_dist = dist;
                    min_id = bn->id;
                }
            }
            bn = bn->next;
        }
        ij = neighbors[ci];
        ci++;
        if(ci == num_neighbors) break;
    } while(1);
    return min_id;
}

struct solution nearest_neighbor_solution(const int* ids, int n, int starting_pos)
{
    //delete_bucket_grid_2d(bucket_grid);
    //int bucket_grid_len = GRID_LEN;
    //struct bucket_grid_2d bucket_grid = init_bucket_grid_2d(bucket_grid_len);
    for(int j = 0; j < n; j++)
    {
        int i = ids[j];
        struct pos ij = get_pos_from_xy(cities[i].x, cities[i].y);
        insert_at_bucket_grid_2d(bucket_grid, ij.i, ij.j, cities[i].id);
    }
    struct solution ret;
    ret.n = n;
    ret.ids = malloc(n * sizeof(int));
    ret.ids[0] = starting_pos;
    ret.ids[n-1] = starting_pos;
    int max_id = 0;
    for(int i = 0; i < n; i++)
    {
        if(ids[i] > max_id) max_id = ids[i];
    }
    char* visited = calloc(max_id + 1, 1);
    visited[starting_pos] = 1;
    for(int i = 1; i < n - 1; i++)
    {
       int next_id = get_nearest_neighbor(ret.ids[i - 1], visited, bucket_grid, 0);//i % 10 == 0);
       assert(next_id != starting_pos);
       visited[next_id] = 1;
       ret.ids[i] = next_id;
       if(i % 10000 == 0) printf("%lf\n", (double)i / n);
    }
    printf("Last Edge Cost: %lf\n", dist_between_id(ret.ids[n-1], ret.ids[n-2]));
    return ret;
}

struct solution nearest_neighbor_dual_approach_solution(const int* ids, int n)
{
    //delete_bucket_grid_2d(bucket_grid);
    //int bucket_grid_len = GRID_LEN;
    //struct bucket_grid_2d bucket_grid = init_bucket_grid_2d(bucket_grid_len);
    for(int j = 0; j < n; j++)
    {
        int i = ids[j];
        struct pos ij = get_pos_from_xy(cities[i].x, cities[i].y);
        insert_at_bucket_grid_2d(bucket_grid, ij.i, ij.j, cities[i].id);
    }
    struct solution ret;
    ret.n = n;
    ret.ids = malloc(n * sizeof(int));
    ret.ids[0] = 0;
    ret.ids[n-1] = 0;
    int max_id = 0;
    for(int i = 0; i < n; i++)
    {
        if(ids[i] > max_id) max_id = ids[i];
    }
    char* visited = calloc(max_id + 1, 1);
    visited[0] = 1;
    int i, j;
    i = 1;
    j = n - 2;
    while(i <= j)
    {
       int next_id = get_nearest_neighbor(ret.ids[i - 1], visited, bucket_grid, i % 10 == 0);
       assert(next_id != 0);
       visited[next_id] = 1;
       ret.ids[i] = next_id;
       if(i == j) break;
       next_id = get_nearest_neighbor(ret.ids[j + 1], visited, bucket_grid, (j + 1) % 10 == 0);
       assert(next_id != 0);
       visited[next_id] = 1;
       ret.ids[j] = next_id;

       if(i % 10000 == 0) printf("%lf\n", (i * 2.0)/ n);
       i++;
       j--;
    }
    printf("Last Edge Cost: %lf\n", dist_between_id(ret.ids[n-1], ret.ids[n-2]));
    return ret;
}


int mutate_h3_solution(struct solution sol, int extent);
double eval_solution(struct solution sol);
struct solution generate_solution(const int* ids, int n)
{
    struct solution ret;
    //ret.n = n;
    //ret.ids = malloc(n * sizeof(int));
    //memcpy(ret.ids, ids, n * sizeof(int));
    //random_permutattion(ret.ids + 1, n - 2);
    //int max_step = 10000;
    //while(max_step--) mutate_h3_solution(ret, 100);
    ret = nearest_neighbor_solution(ids, n, STARTING_POS);
    double cost = eval_solution(ret);
    printf("Initial solution cost: %lf\n", cost);
    return ret;
}

void mutate_solution(struct solution sol)
{
    int n = sol.n - 2;
    int pos1 = 1 + rand() % n;
    int pos2 = 1 + rand() % n;
    int temp = sol.ids[pos1];
    sol.ids[pos1] = sol.ids[pos2];
    sol.ids[pos2] = temp;
}

int mutate_h1_solution(struct solution sol)
{
    int ret = 0;
    int i = rand() % (sol.n - 3);
    int k = 3 + rand() % 500;
    if(i + k > sol.n) k -= sol.n - i;
    int a = sol.ids[i];
    int b = sol.ids[i + 1];
    int c = sol.ids[i + k - 1];
    int d = sol.ids[i + k];
    if(dist_between_id(a, c) + dist_between_id(b, d) < dist_between_id(a, b) + dist_between_id(c, d))
    {
        int j = i + k - 1;
        i++;
        while(i <= j)
        {
            int temp = sol.ids[i];
            sol.ids[i] = sol.ids[j];
            sol.ids[j] = temp;
            i++;
            j--;
        }
        ret = 1;
    }
    return ret;
}

void mutate_h2_solution(struct solution sol)
{
    if(mutate_h1_solution(sol));
    else mutate_solution(sol);
}

int mutate_h3_solution(struct solution sol, int extent)
{
    int ret = 0;
    int n = sol.n;
    if(rand() / RAND_MAX > 0.5)
    {
        for(int i = 0; i < n - 2; i++)
        {
            double prev_dist = dist_between_id(sol.ids[i], sol.ids[i + 1]);
            for(int j = i + 2; j < n - 1 && j < i + extent; j++)
            {
                if(dist_between_id(sol.ids[i], sol.ids[j]) < prev_dist)
                {
                    int temp = sol.ids[i + 1];
                    sol.ids[i + 1] = sol.ids[j];
                    sol.ids[j] = temp;
                    ret = 1;
                    break;
                }
            }
        }
    }
    else
    {
        for(int i = n - 1; i > 1; i--)
        {
            double prev_dist = dist_between_id(sol.ids[i], sol.ids[i - 1]);
            for(int j = i - 1; j > 0 && j > i - extent; j--)
            {
                if(dist_between_id(sol.ids[i], sol.ids[j]) < prev_dist)
                {
                    int temp = sol.ids[i - 1];
                    sol.ids[i - 1] = sol.ids[j];
                    sol.ids[j] = temp;
                    ret = 1;
                    break;
                }
            }
        }
    }

    return ret;
}

void mutate_h4_solution(struct solution sol)
{

}

void multiple_mutate_solution(struct solution sol, int num)
{
    for(int i = 0; i < num; i++)
    {
        mutate_h1_solution(sol);
    }
}

struct solution copy_solution(struct solution sol)
{
    struct solution ret;
    ret.n = sol.n;
    ret.ids = malloc(sol.n * sizeof(int));
    memcpy(ret.ids, sol.ids, sol.n * sizeof(int));
    return ret;
}

void inplace_copy_solution(struct solution dest, struct solution src)
{
    memcpy(dest.ids, src.ids, src.n * sizeof(int));
}

void delete_solution(struct solution sol)
{
    free(sol.ids);
}

double eval_solution(struct solution sol)
{
    double cost = 0.0;
    int i;
    #pragma omp parallel for reduction (+:cost)
    for(i = 0; i < sol.n - 1; i++)
    {
        double d = dist_between_id(sol.ids[i], sol.ids[i + 1]);
        if(i + 1 % 10 == 0 && !is_prime[sol.ids[i + 1]]) d *= 1.1;
        cost += d;
    }
    return cost;
}

void assert_solution(struct solution sol)
{
    /*
    assert(sol.n > 0);
    assert(sol.ids[0] == 0);
    assert(sol.ids[sol.n - 1] == 0);
    */
}

struct solution steepest_ascent_hc(int neighbor_count, long long num_iter, int multi_mut, double gamma, int patience, struct solution template_solution)
{
    int best_sol_index = 0;
    //struct solution best_solution = generate_solution(template_solution.ids, template_solution.n);
    struct solution best_solution = copy_solution(template_solution);
    double best_cost = eval_solution(best_solution);
    printf("Initial cost: %lf\n", best_cost);
    struct solution* sols;
    sols = calloc(neighbor_count, sizeof(struct solution));
    for(int j = 0; j < neighbor_count; j++)
    {
        sols[j] = copy_solution(best_solution);
    }

    int no_improvement = 0;
    //int mut3_extent = 10;
    //int should_use_h3 = 1;
    for(long long i = 0; no_improvement < patience && i < num_iter; i++)
    {
        // make copies of best solution
        for(int j = 0; j < neighbor_count; j++)
        {
            inplace_copy_solution(sols[j], best_solution);
        }
        //puts("made copies of best sol");
        // mutate each solution
        for(int j = 0; j < neighbor_count; j++)
        {
            //mutate_solution(sols[j]);
            multiple_mutate_solution(sols[j], (int)ceil(multi_mut * exp(gamma * i)));
            //mutate_h2_solution(sols[j]);
            //mutate_h1_solution(sols[j]);
            //if(should_use_h3) should_use_h3 = mutate_h3_solution(sols[j], mut3_extent);
            //else mutate_h2_solution(sols[j]);
        }
        // find best solution
        best_sol_index = -1;
        for(int j = 0; j < neighbor_count; j++)
        {
            double cost = eval_solution(sols[j]);
            if(cost < best_cost)
            {
                best_sol_index = j;
                best_cost = cost;
            }
        }
        if(best_sol_index == -1)
        {
            //puts("did not find better solution");
            // got only worse solutions
            // copy previous best solution to position 0
            best_sol_index = 0;
            inplace_copy_solution(sols[0], best_solution);
            no_improvement++;
        }
        else
        {
            //puts("found better solution");
            // copy best soluion to best_solution
            inplace_copy_solution(best_solution, sols[best_sol_index]);
            no_improvement = 0;
        }
        if(i % 100 == 0) printf("Iter: %lld Best Cost: %lf\n", i, best_cost);
    }
    for(int j = 0; j < neighbor_count; j++)
        delete_solution(sols[j]);
    free(sols);
    return best_solution;
}

int main(int argc, char* argv[])
{
    srand(time(NULL));
    load_cities(argv[1]);
    //int starting_pos;
    //sscanf(argv[3], "%d", &starting_pos);
    //printf("%d\n", starting_pos);
    int* primes;
    int prime_count = get_all_primes(N, &primes);
    //int n = prime_count + 2;
    int n = N;
    printf("n = %d\n", n);
    //int bucket_grid_len = GRID_LEN;
    //bucket_grid = init_bucket_grid_2d(bucket_grid_len);
    //for(int i = 0; i < n; i++)
    //{
    //    struct pos ij = get_pos_from_xy(cities[i].x, cities[i].y);
    //    insert_at_bucket_grid_2d(bucket_grid, ij.i, ij.j, cities[i].id);
    //}
    //puts("Bucket grid creation done");

    struct solution template_solution = load_solution(argv[3], n);
    //template_solution.n = n;
    //template_solution.ids = malloc(n * sizeof(int));
    //template_solution.ids[0] = STARTING_POS;
    //template_solution.ids[n - 1] = STARTING_POS;
    //for(int i = 1; i < n - 1; i++)
    //{
        //template_solution.ids[i] = primes[i - 1];
    //    template_solution.ids[i] = i;
    //}
    //struct solution result = nearest_neighbor_solution(template_solution.ids, n, starting_pos);
    //struct solution result = nearest_neighbor_dual_approach_solution(template_solution.ids, n);
    struct solution result = steepest_ascent_hc(20, 10000LL, 100000, -0.00001, 500, template_solution);
    //puts("Done with nearest neighbor");
    double cost = eval_solution(result);
    printf("Cost: %lf\n", cost);
    int pos_of_0 = 0;
    //for(int i = 0; i < result.n; i++)
    //{
    //    if(result.ids[i] == 0)
    //    {
    //        pos_of_0 = i;
    //        break;
    //    }
    //}
    FILE* f = fopen(argv[2], "w");
    fprintf(f, "Path\n");
    for(int i = pos_of_0; i < result.n; i++)
    {
        fprintf(f, "%d\n", result.ids[i]);
    }
    for(int i = 1; i <= pos_of_0; i++)
    {
        fprintf(f, "%d\n", result.ids[i]);
    }
    fclose(f);
    delete_solution(template_solution);
    delete_solution(result);
    //delete_bucket_grid_2d(bucket_grid);
    free(primes);
    return 0;
}
