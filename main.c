#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>

#define N 197770 

int get_all_primes(int below, int** out)
{
    int* ret = calloc(below, sizeof(int));
    // ret[p] == 0 means p is a prime
    if (ret == NULL)
        return -1;
    ret[0] = 1;
    ret[1] = 1;
    for(int i = 2; i < below; i++)
    {
        if(!ret[i])
        {
            for(int q = 2 * i; q < below; q += i)
                ret[q] = 1;
        }
    }
    int current_prime_index = 0;
    for(int i = 0; i < below; i++)
    {
        if(!ret[i])
            ret[current_prime_index++] = i;
    }
    ret = realloc(ret, current_prime_index * sizeof(int));
    *out = ret;
    return current_prime_index;
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


struct solution
{
    int n;
    int* ids;
};

int mutate_h3_solution(struct solution sol, int extent);
struct solution generate_solution(const int* ids, int n)
{
    struct solution ret;
    ret.n = n;
    ret.ids = malloc(n * sizeof(int));
    memcpy(ret.ids, ids, n * sizeof(int));
    random_permutattion(ret.ids + 1, n - 2);
    int max_step = 10000;
    while(max_step--) mutate_h3_solution(ret, 100);
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
    int max_attempt = sol.n;
    int ret = 0;
    while(max_attempt--)
    {
        int i = rand() % (sol.n - 3);
        int a = sol.ids[i];
        int b = sol.ids[i + 1];
        int c = sol.ids[i + 2];
        int d = sol.ids[i + 3];
        if(dist_between_id(a, c) + dist_between_id(b, d) < dist_between_id(a, b) + dist_between_id(c, d))
        {
            sol.ids[i + 1] = c;
            sol.ids[i + 2] = b;
            ret = 1;
        }
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

void multiple_mutate_solution(struct solution sol, int num)
{
    for(int i = 0; i < num; i++)
    {
        mutate_solution(sol);
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
    struct solution best_solution = generate_solution(template_solution.ids, template_solution.n);
    double best_cost = eval_solution(best_solution);
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
        if(i % 10 == 0) printf("Iter: %lld Best Cost: %lf\n", i, best_cost);
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
    int* primes;
    int prime_count = get_all_primes(N, &primes);
    printf("%d\n", prime_count);
    //int n = prime_count + 2;
    int n = N;
    struct solution template_solution;
    template_solution.n = n;
    template_solution.ids = malloc(n * sizeof(int));
    template_solution.ids[0] = 0;
    template_solution.ids[n - 1] = 0;
    for(int i = 1; i < n - 1; i++)
    {
        //template_solution.ids[i] = primes[i - 1];
        template_solution.ids[i] = i;
    }
    struct solution result = steepest_ascent_hc(1000, 1000LL, 1000, -0.01, 20, template_solution);
    delete_solution(template_solution);
    delete_solution(result);
    free(primes);
    return 0;
}
