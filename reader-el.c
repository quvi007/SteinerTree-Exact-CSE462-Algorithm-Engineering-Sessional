#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#include<assert.h>
#include<sys/utsname.h>
#include<math.h>
#include<ctype.h>

/************************************************************* Configuration. */
#define BIN_HEAP
#define TRACK_OPTIMAL

#ifdef TRACK_RESOURCES
#include<omp.h>
#define TRACK_MEMORY
#define TRACK_BANDWIDTH
#endif

#define MAX_K 32

typedef int index_t; // default to 32-bit indexing

/********************************************************** Global constants. */

#define MAX_DISTANCE ((index_t)0x7FFFFFFF)
#define MATH_INF ((index_t)0x7FFFFFFF)
#define UNDEFINED -1

/************************************************************* Common macros. */

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

// Linked list navigation macros. 

#define pnlinknext(to,el) { (el)->next = (to)->next; (el)->prev = (to); (to)->next->prev = (el); (to)->next = (el); }
#define pnlinkprev(to,el) { (el)->prev = (to)->prev; (el)->next = (to); (to)->prev->next = (el); (to)->prev = (el); }
#define pnunlink(el) { (el)->next->prev = (el)->prev; (el)->prev->next = (el)->next; }
#define pnrelink(el) { (el)->next->prev = (el); (el)->prev->next = (el); }


/*********************************************************** Error reporting. */

#define ERROR(...) error(__FILE__,__LINE__,__func__,__VA_ARGS__);

static void error(const char *fn, int line, const char *func, 
                  const char *format, ...) 
{
    va_list args;
    va_start(args, format);
    fprintf(stderr, 
            "ERROR [file = %s, line = %d] "
            "%s: ",
            fn,
            (index_t)line,
            func);
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);
    abort();    
}

/********************************************************* Get the host name. */

#ifdef TRACK_RESOURCES
#define MAX_HOSTNAME 256

const char *sysdep_hostname(void)
{
    static char hn[MAX_HOSTNAME];

    struct utsname undata;
    uname(&undata);
    strcpy(hn, undata.nodename);
    return hn;
}
#endif

#ifdef TRACK_MEMORY
#define MALLOC(x) malloc_wrapper(x)
#define CALLOC(x, y) calloc_wrapper((x), (y))
#define FREE(x) free_wrapper(x)

#else

#define MALLOC(x) malloc((x))
#define CALLOC(x, y) calloc((x), (y))
#define FREE(x) free((x))

#endif

#ifdef TRACK_MEMORY
index_t malloc_balance = 0;

struct malloc_track_struct
{
    void *p;
    size_t size;
    struct malloc_track_struct *prev;
    struct malloc_track_struct *next;
};

typedef struct malloc_track_struct malloc_track_t;

malloc_track_t malloc_track_root;
size_t malloc_total = 0;

#define MEMTRACK_STACK_CAPACITY 512
size_t memtrack_stack[MEMTRACK_STACK_CAPACITY];
index_t memtrack_stack_top = -1;

void *malloc_wrapper(size_t size)
{
    if(malloc_balance == 0) {
        malloc_track_root.prev = &malloc_track_root;
        malloc_track_root.next = &malloc_track_root;
    }
    void *p = malloc(size);
    if(p == NULL)
        ERROR("malloc fails");
    malloc_balance++;

    malloc_track_t *t = (malloc_track_t *) malloc(sizeof(malloc_track_t));
    t->p = p;
    t->size = size;
    pnlinkprev(&malloc_track_root, t);
    malloc_total += size;
    for(index_t i = 0; i <= memtrack_stack_top; i++)
        if(memtrack_stack[i] < malloc_total)
            memtrack_stack[i] = malloc_total;
    return p;
}

void *calloc_wrapper(size_t n, size_t size)
{
    if(malloc_balance == 0) {
        malloc_track_root.prev = &malloc_track_root;
        malloc_track_root.next = &malloc_track_root;
    }
    void *p = calloc(n, size);
    if(p == NULL)
        ERROR("malloc fails");
    malloc_balance++;

    malloc_track_t *t = (malloc_track_t *) malloc(sizeof(malloc_track_t));
    t->p = p;
    t->size = (n*size);
    pnlinkprev(&malloc_track_root, t);
    malloc_total += (n*size);
    for(index_t i = 0; i <= memtrack_stack_top; i++)
        if(memtrack_stack[i] < malloc_total)
            memtrack_stack[i] = malloc_total;
    return p;
}

void free_wrapper(void *p)
{
    malloc_track_t *t = malloc_track_root.next;
    for(;
        t != &malloc_track_root;
        t = t->next) {
        if(t->p == p)
            break;
    }
    if(t == &malloc_track_root)
        ERROR("FREE issued on a non-tracked pointer %p", p);
    malloc_total -= t->size;
    pnunlink(t);
    free(t);

    free(p);
    malloc_balance--;
}

index_t *alloc_idxtab(index_t n)
{
    index_t *t = (index_t *) MALLOC(sizeof(index_t)*n);
    return t;
}

void push_memtrack(void)
{
    assert(memtrack_stack_top + 1 < MEMTRACK_STACK_CAPACITY);
    memtrack_stack[++memtrack_stack_top] = malloc_total;
}

size_t pop_memtrack(void)
{
    assert(memtrack_stack_top >= 0);
    return memtrack_stack[memtrack_stack_top--];
}

size_t current_mem(void)
{
    return malloc_total;
}

double inGiB(size_t s)
{
    return (double) s / (1 << 30);
}

void print_current_mem(void)
{
    fprintf(stdout, "{curr: %.2lfGiB}", inGiB(current_mem()));
    fflush(stdout);
}

void print_pop_memtrack(void)
{
    fprintf(stdout, "{peak: %.2lfGiB}", inGiB(pop_memtrack()));
    fflush(stdout);
}

void inc_malloc_total(size_t size)
{
    malloc_total += size;
}

void dec_malloc_total(size_t size)
{
    malloc_total -= size;
}
#endif

/******************************************************** Timing subroutines. */

#ifdef TRACK_RESOURCES
#define TIME_STACK_CAPACITY 256
double start_stack[TIME_STACK_CAPACITY];
index_t start_stack_top = -1;

void push_time(void)
{
    assert(start_stack_top + 1 < TIME_STACK_CAPACITY);
    start_stack[++start_stack_top] = omp_get_wtime();
}

double pop_time(void)
{
    double wstop = omp_get_wtime();
    assert(start_stack_top >= 0);
    double wstart = start_stack[start_stack_top--];
    return (double) (1000.0*(wstop-wstart));
}
#endif

/********************************************************* Utility functions. */

/******************************************************* String manipulation. */
char* strlower( char *s)
{
   char* t;

   for(t = s; *s != '\0'; s++)
      *s = (char)tolower(*s);

   return(t);
}

/**************************************************************  prefix sum. */
index_t prefixsum(index_t n, index_t *a, index_t k)
{
    index_t run = 0;
    for(index_t u = 0; u < n; u++) {
        index_t tv = a[u];
        a[u] = run;
        run += tv + k;
    }
    return run;
}


/*************************************************** Graph build subroutines. */

typedef struct graph
{
    index_t n;
    index_t m;
    index_t k;
    index_t num_edges;
    index_t num_terminals;
    index_t edge_capacity;
    index_t cost;
    index_t *edges;
    index_t *terminals;
} graph_t;

static index_t *enlarge(index_t m, index_t m_was, index_t *was)
{
    assert(m >= 0 && m_was >= 0);

    index_t *a = (index_t *) MALLOC(sizeof(index_t)*m);
    index_t i;
    if(was != (void *) 0) { 
        for(i = 0; i < m_was; i++) {
            a[i] = was[i];
        }
        FREE(was);
    }    
    return a;
}

graph_t *graph_alloc()
{
    graph_t *g = (graph_t *) MALLOC(sizeof(graph_t));
    g->n = 0; 
    g->m = 0; 
    g->k = 0;
    g->num_edges      = 0;
    g->num_terminals  = 0;
    g->edge_capacity  = 100;
    g->edges          = enlarge(3*g->edge_capacity, 0, (void *) 0);
    g->terminals      = NULL;
    g->cost           = -1;
    
    return g;
}

void graph_free(graph_t *g)
{
    if(g->edges != NULL)
        FREE(g->edges);
    if(g->terminals != NULL)
        FREE(g->terminals);
    FREE(g);
}

void graph_add_edge(graph_t *g, index_t u, index_t v, index_t w)
{
    assert(u >= 0 && v >= 0 && u < g->n && v < g->n);

    if(g->num_edges == g->edge_capacity)
    {
        g->edges = enlarge(6*g->edge_capacity, 3*g->edge_capacity, g->edges);
        g->edge_capacity *= 2;
    }

    assert(g->num_edges < g->edge_capacity);

    index_t *e = g->edges + 3*g->num_edges;
    g->num_edges++;
    e[0] = u;
    e[1] = v;
    e[2] = w;
}

void graph_add_terminal(graph_t *g, index_t u)
{
    if(g->terminals == NULL)
        ERROR("section terminals not initialised");

    assert(u >= 0 && u < g->n);
    index_t *t = g->terminals + g->num_terminals;
    g->num_terminals++;

    assert(g->num_terminals <= g->k);
    t[0] = u;
}

#define MAX_LINE_SIZE 1024

graph_t * graph_load(FILE *in)
{
#ifdef TRACK_RESOURCES
    push_time();
    push_memtrack();
#endif

    char buf[MAX_LINE_SIZE];
    char in_line[MAX_LINE_SIZE];
    index_t n = 0;
    index_t m = 0;
    index_t k = 0;
    index_t u, v, w;
    index_t cost = -1;
    graph_t *g = graph_alloc();

    while(fgets(in_line, MAX_LINE_SIZE, in) != NULL)
    {
        char *line = strlower(in_line);
        int c = line[0];
        char *tok;
        strcpy(buf, line);
        tok = strtok(buf, " ");
        switch(c) {
        case 'c':
            if(!strcmp(tok, "cost")) {
                sscanf(line, "cost %d", &cost);
                g->cost = cost;
            }
            break;
        case 'e':
            if(!strcmp(tok, "edges")) {
                sscanf(line, "edges %d", &m);
                g->m = m;
                break;
            } 
            else if(!strcmp(tok, "e")) {
                sscanf(line, "e %d %d %d", &u, &v, &w);
                graph_add_edge(g, u-1, v-1, w);
                break;
            }
            break;
        case 'n':
            if(!strcmp(tok, "nodes")) {
                sscanf(line, "nodes %d", &n);
                g->n = n;
            }
            break;
        case 't':
            if(!strcmp(tok, "terminals")) {
                sscanf(line, "terminals %d", &k);
                g->k = k;
                g->terminals = (index_t *) MALLOC(k*sizeof(index_t));
                break;
            }
            if(!strcmp(tok, "t")) {
                sscanf(line, "t %d", &u);
                graph_add_terminal(g, u-1);
                break;
            }
            break;
        default:
            break;
        }
    }

    assert(g->n != 0);
    assert(g->m == g->num_edges && g->m != 0);
    assert(g->k == g->num_terminals && g->k != 0);

#ifdef TRACK_RESOURCES
    double time = pop_time();
    fprintf(stdout, "input: n = %d, m = %d, k = %d, cost = %d [%.2lfms] ",
                    g->n, g->m, g->k, g->cost, time);
    print_pop_memtrack();
    fprintf(stdout, " ");
    print_current_mem();
    fprintf(stdout, "\n");
    fprintf(stdout, "terminals:");
    for(index_t i=0; i<g->k; i++) 
        fprintf(stdout, " %d", g->terminals[i]+1);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    return g;
}




/******************************************************** Root query builder. */

/***************************************************************** Indexing. */
// subset major index
#define FV_INDEX(v, n, k, X) (((index_t)((X)-1) * (n)) + (v))
#define BV_INDEX(v, n, k, X) (((index_t)((X)-1) * (2*(n))) + (2*(v)))


typedef struct steinerq
{
    index_t     n;
    index_t     m;
    index_t     k;
    index_t     *kk;
    index_t     *pos;
    index_t     *adj;
}steinerq_t;

steinerq_t *root_build(graph_t *g)
{
#ifdef TRACK_RESOURCES
    fprintf(stdout, "root build: ");
    push_time();
    push_memtrack();
#endif
    index_t n = g->n;
    index_t m = g->m;
    index_t k = g->k;
    index_t *kk = (index_t *) MALLOC(k*sizeof(index_t));
    index_t *pos = (index_t *) MALLOC((n+1)*sizeof(index_t));
    index_t *adj = (index_t *) MALLOC(((n+1)+(4*m)+(2*n))*sizeof(index_t));

    steinerq_t *root = (steinerq_t *) MALLOC(sizeof(steinerq_t));
    root->n = n;
    root->m = m;
    root->k = k;
    root->kk  = kk;
    root->pos = pos;
    root->adj = adj;

#ifdef TRACK_RESOURCES
    push_time();
#endif
    for(index_t u = 0; u < n; u++)
        pos[u] = 0;
#ifdef TRACK_RESOURCES
    double time = pop_time();
    fprintf(stdout, "[zero: %.2lfms] ", time);
    fflush(stdout);
    push_time();
#endif

    index_t *e = g->edges;
    for(index_t j = 0; j < 3*m; j+=3)
    {
        pos[e[j]]+=2;
        pos[e[j+1]]+=2;
    }

    pos[n] = (2*n);
    index_t run = prefixsum(n+1, pos, 1);
    assert(run == ((n+1)+(4*m)+(2*n)));

#ifdef TRACK_RESOURCES
    time = pop_time();
    fprintf(stdout, "[pos: %.2lfms] ", time);
    fflush(stdout);
    push_time();
#endif

    for(index_t u = 0; u < n; u++)
        adj[pos[u]] = 0;

    for(index_t j = 0; j < 3*m; j+=3)
    {
        index_t u    = e[j+0];
        index_t v    = e[j+1];
        index_t w    = e[j+2];
        index_t pu   = pos[u];
        index_t pv   = pos[v];
        index_t i_pu = pu + 1 + (2*adj[pu]);
        index_t i_pv = pv + 1 + (2*adj[pv]);

        adj[i_pv]     = u;
        adj[i_pv + 1] = w;
        adj[pv]++;
        adj[i_pu]     = v;
        adj[i_pu + 1] = w;
        adj[pu]++; 
    }

    index_t u = n;
    adj[pos[u]] = 0;
    index_t pu = pos[u];
    for(index_t v = 0; v < n; v++)
    {
        index_t i_pu  = pu + 1 + (2*adj[pu]);
        adj[i_pu]     = v;
        adj[i_pu + 1] = MATH_INF;
        adj[pu]++;
    }

    index_t *tt = g->terminals;
#ifdef TRACK_RESOURCES
    pop_time();
    fprintf(stdout, "[adj: %.2lfms] ", time);
    fflush(stdout);
	push_time();
#endif

    for(index_t u = 0; u < k; u++)
        kk[u] = tt[u];

#ifdef TRACK_RESOURCES
    time = pop_time();
    fprintf(stdout, "[term: %.2lfms] ", time);
    fflush(stdout);

    time = pop_time();
    fprintf(stdout, "done. [%.2lfms] ", time);
    print_pop_memtrack();
    fprintf(stdout, " ");
    print_current_mem();
    fprintf(stdout, "\n");
    fflush(stdout);
#endif

    return root;
}

void steinerq_free(steinerq_t *root)
{
    if(root->pos != NULL)
        FREE(root->pos);
    if(root->adj != NULL)
        FREE(root->adj);
    if(root->kk != NULL)
        FREE(root->kk);
    FREE(root);
}

/********************************************************* Debug routines. */

#ifdef DEBUG
void print_nbits(index_t n, index_t N)
{            
    index_t bits[64];
    for(index_t i=0; i<64; i++)
        bits[i] = (N>>i)&0x01;
    for(index_t i=n-1; i >= 0; i--)
        fprintf(stdout, "%d", bits[i]);
    fflush(stdout);
}

void print_bits(index_t n)
{
    fprintf(stdout, "n: %d bits: ", n); 
    index_t size = sizeof(index_t)*8;

    for(index_t i = 0; i<size; i++)
    {   
        index_t mask = ((index_t)(0x01) << (size-1-i));
        fprintf(stdout, "%d", (index_t)((n&mask) ? 1 : 0));
    }
    fprintf(stdout, "\n");
    fflush(stdout);
}

void print_array(index_t n, index_t *a)
{
    fprintf(stdout, "n: %d a:", n);

    for(index_t i = 0; i < n; i++)
    {
        fprintf(stdout, " %d", a[i]+1);
    }
    fprintf(stdout, "\n");
    fflush(stdout);
}

void print_adj(index_t u, index_t *pos, index_t *adj)
{
    index_t p = pos[u];
    index_t nu = adj[p];
    index_t *adj_u = adj + p + 1;
    fprintf(stdout, "adjacency list (%d) : ", u);
    for(index_t i = 0; i < nu; i++)
        fprintf(stdout, " %d %d|", adj_u[2*i]+1, adj_u[2*i+1]);
    fprintf(stdout, "\n");
}

void print_dist(index_t n, index_t *d)
{
    fprintf(stdout, "Shortest distance: \n");
    for(index_t u = 0; u < n; u++)
        fprintf(stdout, "%d: %d\n", u+1, d[u]);
    fflush(stdout);
}


void print_dist_matrix(index_t *d_N, index_t n)
{
    fprintf(stdout, "distance matrix: \n");
    for(index_t u = 0; u < n; u++)
    {
        for(index_t v = 0; v < n; v ++)
        {
            fprintf(stdout, " %d", d_N[u*n+v]);
        }
        fprintf(stdout, "\n");
    }
    fflush(stdout);
}

void print_graph(graph_t *g)
{
    fprintf(stdout, "graph_t: \n");
    fprintf(stdout, "n: %d\n", g->n);
    fprintf(stdout, "m: %d\n", g->m);
    fprintf(stdout, "k: %d\n", g->k);
    fprintf(stdout, "num edges: %d\n", g->num_edges);
    fprintf(stdout, "num terminals: %d\n", g->num_terminals);
    fprintf(stdout, "edge capacity: %d\n", g->edge_capacity);
    fprintf(stdout, "edges: \n");
    for(index_t i = 0; i < g->num_edges; i++) 
    {    
        index_t *e = g->edges + (3*i);
        index_t u = e[0];
        index_t v = e[1];
        index_t w = e[2];
        fprintf(stdout, "E %d %d %d\n", u+1, v+1, w);
    }

    fprintf(stdout, "terminals: \n");
    for(index_t i = 0; i < g->num_terminals; i++) 
        fprintf(stdout, "T %d\n", g->terminals[i]+1);

    fflush(stdout);
}

void print_steinerq(steinerq_t *root)
{
    fprintf(stdout, "steinerq_t: \n");
    fprintf(stdout, "n: %d\n", root->n);
    fprintf(stdout, "m: %d\n", root->m);
    fprintf(stdout, "k: %d\n", root->k);
    index_t *pos = root->pos;
    index_t *adj = root->adj;
    fprintf(stdout, "pos:");
    for(index_t i = 0; i < root->n; i++)
        fprintf(stdout, " %d", pos[i]);
    fprintf(stdout, "\nadj:\n");
    index_t n = root->n + 1;
    for(index_t u = 0; u < n; u++)
    {
        index_t pu = pos[u];
        index_t adj_u = adj[pu];
        fprintf(stdout, "node: %d edges: %d|", u+1, adj_u);
        for(index_t i = 0; i < adj_u; i++)
        {
            fprintf(stdout, " %d %d|", 
                            adj[pu + 1 + (2*i)]+1, 
                            adj[pu + 1 + (2*i+1)]);
        }
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
    fflush(stdout);
}

void print_f_v(index_t n, index_t k, index_t *f_v)
{
    fprintf(stdout, "f_v: \n");
    for(index_t X = 0; X < (1<<(k-1)); X++) 
    {    
        fprintf(stdout, "X->");
        print_nbits(k, X+1);
        fprintf(stdout, "\n");
     
        index_t *f_X = f_v + (X*n);
        for(index_t v = 0; v < n; v++) 
        {
            fprintf(stdout, "%d", v+1);
            if(f_X[v] == MATH_INF)
                fprintf(stdout, " MATH_INF\n");
            else
                fprintf(stdout, " %d\n", f_X[v]);
        }
    }
    fflush(stdout);
}

void print_b_v(index_t n, index_t k, index_t *b_v)
{
    fprintf(stdout, "b_v: \n");
    for(index_t X = 0; X < (1<<(k-1)); X++)
    {
        fprintf(stdout, "X->\n");
        print_nbits(k, X+1);
        fprintf(stdout, "\n");

        index_t *b_X = b_v + ((index_t)(2)*X*n);

        for(index_t v = 0; v < n; v++)
        {
            index_t u = b_X[2*v];
            index_t Xd = b_X[2*v + 1];

            fprintf(stdout, "v: %d u: %d v: ",
                    v==-1?-1:v+1, u==-1?-1:u+1);
            print_nbits(k, Xd);
            fprintf(stdout, "\n");
        }
    }
	fflush(stdout);
}
#endif 

/******************************************************* Heap implementaions. */

/************************************************ Binary heap implementation. */
#ifdef BIN_HEAP
typedef struct bheap_item
{
    index_t item;
    index_t key; 
} bheap_item_t;

typedef struct bheap
{
    index_t max_n;
    index_t n;       // size of binary heap
    bheap_item_t *a; // stores (distance, vertex) pairs of the binary heap
    index_t *p;      // stores the positions of vertices in the binary heap
#ifdef TRACK_BANDWIDTH
    index_t key_comps;
    index_t mem;
#endif
} bheap_t;

bheap_t * bh_alloc(index_t n)
{
    bheap_t *h = (bheap_t *) malloc(sizeof(bheap_t));
    h->max_n = n;
    h->n = 0; 
    h->a = (bheap_item_t *) malloc((n+1)*sizeof(bheap_item_t));
    h->p = (index_t *) malloc(n*sizeof(index_t));
#ifdef TRACK_BANDWIDTH
    h->key_comps  = 0; 
    h->mem = 0;
#endif
#ifdef TRACK_MEMORY
    inc_malloc_total(sizeof(bheap_t) + 
                     ((n+1)*sizeof(bheap_item_t)) +
                     (n*sizeof(index_t)));
#endif
    return h;
}

void bh_free(bheap_t *h)
{
#ifdef TRACK_MEMORY
    index_t n = h->max_n;
    dec_malloc_total(sizeof(bheap_t) + 
                     ((n+1)*sizeof(bheap_item_t)) +
                     (n*sizeof(index_t)));
#endif
    free(h->a);
    free(h->p);
    free(h);
}

/************************************************** Binary heap operations. */

static void bh_siftup(bheap_t *h, index_t p, index_t q)
{
    index_t j = p;
    index_t k = 2 * p;
    bheap_item_t y = h->a[p];
#ifdef TRACK_BANDWIDTH
    index_t mem = 0;
    index_t key_comps = 0;
#endif

    while(k <= q)
    {
        bheap_item_t z = h->a[k];
        if(k < q)
        {
#ifdef TRACK_BANDWIDTH
            mem++;
            key_comps++;
#endif
            if(z.key > h->a[k + 1].key) z = h->a[++k];
        }

#ifdef TRACK_BANDWIDTH
        mem += 2;
        key_comps++;
#endif
        if(y.key <= z.key) break;
        h->a[j] = z;
        h->p[z.item] = j;
        j = k;
        k = 2 * j;
    }

    h->a[j] = y;
    h->p[y.item] = j;

#ifdef TRACK_BANDWIDTH
    h->mem += mem;
    h->key_comps += key_comps;
#endif
}

bheap_item_t bh_min(bheap_t *h)
{
    return (bheap_item_t) h->a[1];
}

static void bh_insert(bheap_t *h, index_t item, index_t key)
{
    index_t i = ++(h->n);
#ifdef TRACK_BANDWIDTH
    index_t mem = 0;
    index_t key_comps = 0;
#endif

    while(i >= 2)
    {
        index_t j = i / 2;
        bheap_item_t y = h->a[j];

#ifdef TRACK_BANDWIDTH
        mem ++;
        key_comps++;
#endif
        if(key >= y.key) break;

        h->a[i] = y;
        h->p[y.item] = i;
        i = j;
    }

    h->a[i].item = item;
    h->a[i].key = key;
    h->p[item] = i;
#ifdef TRACK_BANDWIDTH
    h->mem += mem;
    h->key_comps += key_comps;
#endif
}

static void bh_delete(bheap_t *h, index_t item)
{
    index_t n = --(h->n);
    index_t p = h->p[item];
#ifdef TRACK_BANDWIDTH
    index_t mem = 0;
    index_t key_comps = 0;
#endif

    if(p <= n)
    {
#ifdef TRACK_BANDWIDTH
        key_comps++;
        mem += 2;
#endif
        if(h->a[p].key <= h->a[n + 1].key)
        {
            h->a[p] = h->a[n + 1];
            h->p[h->a[p].item] = p;
            bh_siftup(h, p, n);
        }
        else
        {
            h->n = p - 1;
            bh_insert(h, h->a[n + 1].item, h->a[n+1].key);
            h->n = n;
        }
    }
#ifdef TRACK_BANDWIDTH
    h->mem += mem;
    h->key_comps += key_comps;
#endif
}

static void bh_decrease_key(bheap_t *h, index_t item, index_t new_key)
{
#ifdef TRACK_BANDWIDTH
    index_t mem = 1;
    index_t key_comps = 0;
#endif

    index_t i = h->p[item];
    while(i >= 2)
    {
        index_t j = i / 2;
        bheap_item_t y = h->a[j];

#ifdef TRACK_BANDWIDTH
        mem ++;
        key_comps++;
#endif
        if(new_key >= y.key) break;

        h->a[i] = y;
        h->p[y.item] = i;
        i = j;
    }

    h->a[i].item = item;
    h->a[i].key = new_key;
    h->p[item] = i;
#ifdef TRACK_BANDWIDTH
    h->mem += mem;
    h->key_comps += key_comps;
#endif
}

static index_t bh_delete_min(bheap_t * h)
{    
    bheap_item_t min = (bheap_item_t) h->a[1];
    index_t u = min.item;
    bh_delete((bheap_t *)h, u);
    return u;
}
#endif


/************************************************** Heap wrapper functions. */

#ifdef BIN_HEAP
// allocation
#define heap_alloc(n) bh_alloc((n))
#define heap_free(h) bh_free((bheap_t *)(h));
// heap operations
#define heap_insert(h, v, k) bh_insert((h), (v), (k))
#define heap_delete_min(h) bh_delete_min((h));
#define heap_decrease_key(h, v, k) bh_decrease_key((h), (v), (k));
// fetch structure elements
#define heap_n(h) ((bheap_t *)h)->n;
#define heap_key_comps(h) ((bheap_t *)h)->key_comps;
#define heap_mem(h) (h)->mem;
// heap nodes
#define heap_node_t bheap_item_t
#define heap_t bheap_t
#endif


/************************************************** Dijkstra shortest path*/

void dijkstra(index_t n,
              index_t m, 
              index_t *pos, 
              index_t *adj, 
              index_t s, 
              index_t *d,
              index_t *visit,
              index_t *p
#ifdef TRACK_BANDWIDTH
              ,index_t *heap_ops
#endif
             )
{

    heap_t *h = heap_alloc(n);

    for(index_t v = 0; v < n; v++)
    {
        d[v]     = MAX_DISTANCE; // mem: n
        visit[v] = 0; // mem: n
        //p[v]     = UNDEFINED; // skip intialising
    }
    d[s] = 0;
    p[s] = UNDEFINED;

    for(index_t v = 0; v < n; v++)
        heap_insert(h, v, d[v]);

    //visit and label
    while(h->n > 0)
    {
        index_t u = heap_delete_min(h); 
        visit[u]  = 1;

        index_t pos_u  = pos[u];
        index_t *adj_u = adj + pos_u;
        index_t n_u    = adj_u[0];
        for(index_t i = 1; i <= 2*n_u; i += 2)
        {
            index_t v   = adj_u[i];
            index_t d_v = d[u] + adj_u[i+1];
            if(!visit[v] && d[v] > d_v)
            {
                d[v] = d_v;
                heap_decrease_key(h, v, d_v);
                p[v] = u;
            }
        }
        // mem: 2n+6m
    }

#ifdef TRACK_BANDWIDTH
    *heap_ops = heap_mem(h);
#endif
    heap_free(h);
}

/*************************************************** traceback Steiner tree. */

// build a path using the output of Dijkstra's algorithm
void tracepath(index_t n, index_t s, index_t cost, index_t v, index_t *p)
{
    fprintf(stdout, "VALUE %d\n", cost);
    index_t u = p[v];
    while(u != s)
    {
        //graph_add_edge(g, v, u, 1);
        fprintf(stdout, "%d %d\n", v+1, u+1);
        v = u;
        u = p[v];
    }
    fprintf(stdout, "%d %d\n", v+1, u+1);
}

// handling more general case
void list_solution(graph_t *g)
{
    index_t m  = g->num_edges;
    index_t *e = g->edges;
    index_t i  = 0;

    fprintf(stdout, "solution: [");
    for(i = 0; i < (3*m-3); i+=3)
    {
        index_t u = e[i];
        index_t v = e[i+1];
        //index_t w = e[i+2];
        fprintf(stdout, "\"%d %d\", ", u+1, v+1);
    }
    index_t u = e[i];
    index_t v = e[i+1];
    //index_t w = e[i+2];
    fprintf(stdout, "\"%d %d\"", u+1, v+1);
    fprintf(stdout, "]\n");
    fflush(stdout);
}

void backtrack(index_t n, index_t k, index_t v, 
               index_t X, index_t *kk, index_t *b_v
#ifdef TRACK_RESOURCES
               ,graph_t *g
#endif
              )
{
    if(X == 0 || v == -1)
        return;

    index_t i_X = BV_INDEX(v, n, k, X);
    index_t u   = b_v[i_X];

    if(v != u)
    {
        if(u == -1)
            return;
#ifdef TRACK_RESOURCES
        graph_add_edge(g, v, u, 1);
#else
        fprintf(stdout, "%d %d\n", v+1, u+1);
#endif
        index_t Xd = b_v[i_X+1];
        backtrack(n, k, u, Xd, kk, b_v
#ifdef TRACK_RESOURCES
                  , g
#endif
                  );
    }
    else
    {
        index_t Xd   = b_v[i_X+1];
        index_t X_Xd = (X & ~Xd);
        if(X == Xd)
            return;
        backtrack(n, k, u, Xd, kk, b_v
#ifdef TRACK_RESOURCES
                  , g
#endif
                  );
        backtrack(n, k, u, X_Xd, kk, b_v
#ifdef TRACK_RESOURCES
                  , g
#endif
                  );
    }
}

void build_tree(index_t n, index_t k, index_t cost, 
                index_t *kk, index_t *b_v
#ifdef TRACK_RESOURCES
                , graph_t *g
#endif
               )
{
    index_t c = k-1;
    index_t C = (1<<c)-1;
    index_t q = kk[k-1];

#ifdef TRACK_RESOURCES
    g->n = n;
#else
    fprintf(stdout, "VALUE %d\n", cost);
#endif

    backtrack(n, k, q, C, kk, b_v
#ifdef TRACK_RESOURCES
              ,g
#endif
             );
}

/**************************************************** Erickson Monma Veinott. */

index_t emv_kernel(index_t n, 
                   index_t m, 
                   index_t k, 
                   index_t *pos, 
                   index_t *adj, 
                   index_t *kk, 
                   index_t *d,  
                   index_t *p,
                   index_t *visit,
                   index_t *f_v, 
                   index_t *b_v 
#ifdef TRACK_BANDWIDTH
                   ,index_t *heap_ops
#endif
                   )
{
    // initialisation
    index_t c = k-1;
    index_t q = kk[k-1];
    index_t C = (1<<c)-1;

    for(index_t t = 0; t < k; t++) 
    {    
        dijkstra(n+1, m, pos, adj, kk[t], d, visit, p
#ifdef TRACK_BANDWIDTH
                 ,heap_ops
#endif
                );
        index_t *f_t = f_v + FV_INDEX(0, n, k, 1<<t);
        index_t *b_t = b_v + BV_INDEX(0, n, k, 1<<t);

        for(index_t v = 0; v < n; v++) 
        {
            f_t[v]       = d[v]; 
            b_t[2*v]     = p[v];
            b_t[2*v + 1] = (1<<t);
        }
        // mem: 3*k*n
    }

    // computing a Steiner tree for set K \ {q}
    for(index_t m = 2; m < c; m++) // k-2
    {    
        index_t z = 0;
        // bit twiddling hacks: generating all subsets of size `m`
        for(index_t X = (1<<m)-1;
            X < (1<<c);
            z = X|(X-1), X = (z+1)|(((~z & -~z)-1) >> (__builtin_ctz(X) + 1))) // cCm
        {
            index_t *f_X  = f_v + FV_INDEX(0, n, k, X);
            index_t *b_X  = b_v + BV_INDEX(0, n, k, X);
            // bit twiddling hacks: generating proper subsets of X 
            index_t Xd = 0;
            for(Xd = X & (Xd - X); Xd < (X/2 + 1); Xd = X & (Xd - X)) // 2^{m-1} 
            {
                index_t X_Xd    = (X & ~Xd); // X - X' 
                index_t *f_Xd   = f_v + FV_INDEX(0, n, k, Xd);
                index_t *f_X_Xd = f_v + FV_INDEX(0, n, k, X_Xd);
                for(index_t v = 0; v < n; v++)
                {
                    index_t min_Xd = f_Xd[v] + f_X_Xd[v];
                    if(min_Xd < f_X[v])
                    {
                        f_X[v]       = min_Xd;
                        b_X[2*v]     = v;
                        b_X[2*v + 1] = Xd;
                    }
                    // mem: 3^c*3n
                }
            }

            index_t s      = n;
            index_t ps     = pos[s];
            index_t *adj_s = adj + (ps+1);
            for(index_t u = 0; u < n; u++)
                adj_s[2*u+1] = f_X[u]; // mem: 2^c * n

            for(index_t t = 0; t < k; t++)
            {
                if(!(X & (1<<t)))
                    continue;
                index_t u     = kk[t];
                index_t X_u   = (X & ~(1<<t));
                index_t i_X_u = FV_INDEX(u, n, k, X_u);
                adj_s[2*u+1]  = f_v[i_X_u];
            }

            dijkstra(n+1, m+n, pos, adj, s, d, visit, p
#ifdef TRACK_BANDWIDTH
                     ,heap_ops
#endif
                     );
            for(index_t v = 0; v < n; v++)
            {
                f_X[v]    = d[v];
                index_t u = p[v];
                if(u != s)
                {
                    b_X[2*v]     = u;
                    b_X[2*v + 1] = X;
                }
                // mem: 2^c * 2n 
            }
        }
    }

    // finalise: computing a Steiner tree for {K \ {q}} U {q}
    index_t *f_C = f_v + FV_INDEX(0, n, k, C);
    index_t *b_C = b_v + BV_INDEX(0, n, k, C);
    index_t Xd   = 0; 
    for(Xd = C & (Xd - C); Xd < (C/2 + 1); Xd = C & (Xd - C)) // 2^{k-1}
    {    
        index_t C_Xd    = (C & ~Xd); // C - X'
        index_t *f_Xd   = f_v + FV_INDEX(0, n, k, Xd); 
        index_t *f_C_Xd = f_v + FV_INDEX(0, n, k, C_Xd);
        for(index_t v = 0; v < n ; v++) // n 
        {
            index_t min_Xd = f_Xd[v] + f_C_Xd[v];
            if(min_Xd < f_C[v])
            {
                f_C[v]       = min_Xd;
                b_C[2*v]     = v; 
                b_C[2*v + 1] = Xd;
            }
            // mem: 2^k *3n
        }
    } 
    index_t s      = n;
    index_t ps     = pos[s];
    index_t *adj_s = adj + (ps+1);
    for(index_t u = 0; u < n; u++)
        adj_s[2*u+1] = f_C[u]; // mem: n

    for(index_t t = 0; t < k; t++)
    {
        if(!(C & (1<<t)))
            continue;
        index_t u     = kk[t];
        index_t X_u   = (C & ~(1<<t));
        index_t i_X_u = FV_INDEX(u, n, k, X_u);
        adj_s[2*u+1]  = f_v[i_X_u];
    }

    dijkstra(n+1, m+n, pos, adj, s, d, visit, p
#ifdef TRACK_BANDWIDTH
             ,heap_ops
#endif
             );
    for(index_t v = 0; v < n; v++)
    {
        f_C[v]    = d[v];
        index_t u = p[v];
        if(u != s)
        {
            b_C[2*v]     = u;
            b_C[2*v + 1] = C;
        }
        // mem: 2n
    }

    index_t i_q_C  = FV_INDEX(q, n, k, C);
    return f_v[i_q_C];
}

index_t erickson_monma_veinott(steinerq_t *root)
{
#ifdef TRACK_RESOURCES
    fprintf(stdout, "erickson: ");
    fflush(stdout);
    push_memtrack();
    push_time();
#endif
    index_t n   = root->n;
    index_t m   = root->m;
    index_t k   = root->k;
    index_t *kk = root->kk;
    index_t min_cost = 0;

    if(k > MAX_K)
        ERROR("Maximum supported terminal-set size is '%d'", MAX_K);

    // zero or one terminal
    if(k == 0 || k == 1) 
    {
        fprintf(stdout, "VALUE %d\n", min_cost);
        if(k == 1) {
            fprintf(stdout, "%d\n", kk[0]);
        }
        return min_cost;
    }

    // two terminals
    if(k == 2) 
    {
        index_t u      = kk[0];
        index_t v      = kk[1];
        index_t *d     = (index_t *) MALLOC(n*sizeof(index_t));
        index_t *visit = (index_t *) MALLOC(n*sizeof(index_t));
        index_t *p     = (index_t *) MALLOC(n*sizeof(index_t));
#ifdef TRACK_BANDWIDTH
        index_t heap_ops = 0;
#endif
        dijkstra(n, m, root->pos, root->adj, u, d, visit
                ,p
#ifdef TRACK_BANDWIDTH
                ,&heap_ops
#endif
                );

        // compute bandwidth
        min_cost = d[v];
        tracepath(n, u, min_cost, v, p);

        FREE(d);
        FREE(visit);
        FREE(p);
        return min_cost;
    }

    // handling more general case : more than two terminals
    index_t c = k-1;
    index_t *f_v   = (index_t *) MALLOC(n*(1<<c)*sizeof(index_t));
    index_t *b_v   = (index_t *) MALLOC(2*n*(1<<c)*sizeof(index_t));
    index_t *d     = (index_t *) MALLOC((n+1)*sizeof(index_t));
    index_t *p     = (index_t *) MALLOC((n+1)*sizeof(index_t));
    index_t *visit = (index_t *) MALLOC((n+1)*sizeof(index_t));
#ifdef TRACK_BANDWIDTH
    index_t *heap_ops = (index_t *) MALLOC(sizeof(index_t));
#endif

#ifdef TRACK_RESOURCES
    push_time();
#endif
    // initialisation
    for(index_t i = 0; i < (index_t)(n*(1<<c)); i++)
        f_v[i] = MATH_INF;

    for(index_t i = 0; i < (index_t)(2*n*(1<<c)); i+=2) 
    {
        b_v[i]   = -1;
        b_v[i+1] = 0;
    }

#ifdef TRACK_RESOURCES
    double time = pop_time();
    fprintf(stdout, "[zero: %.2lfms] ", time);
    push_time();
#endif

    // call kernel: do the hard work
    min_cost = emv_kernel(n, m, k, root->pos, root->adj, kk, 
                          d, p, visit, f_v, b_v
#ifdef TRACK_BANDWIDTH
                          ,heap_ops
#endif
                          );
    
#ifdef TRACK_RESOURCES
    time = pop_time();
    index_t total_heap_ops = heap_ops[0];
    index_t mem_graph = ((5*n) + (6*m));
    // mem: ((3^k*n + 2^(k-1)*9n + 3n*(k+1)) + 2^(k-1)*(5n+6m))*sizeof(index_t) 
    //      + total_heap_ops*sizeof(heap_node_t)
    index_t trans_bytes = (((index_t)(pow(3,k)*n)+(index_t)(pow(2,k-1)*9*n)+
                      3*n*(k+1)+(index_t)(pow(2,k-1)*mem_graph))*sizeof(index_t)
                      +total_heap_ops*sizeof(heap_node_t));
    double trans_rate = trans_bytes / (time / 1000.0);
    fprintf(stdout, "[kernel: %.2lfms %.2lfGiB/s] ",
                    time, trans_rate/(1 << 30));
    push_time();
    graph_t *g = graph_alloc();
    g->n = n;
#endif

    // build a Steiner tree
    build_tree(n, k, min_cost, kk, b_v
#ifdef TRACK_RESOURCES
               ,g
#endif
               );

#ifdef TRACK_RESOURCES
    time= pop_time();
    fprintf(stdout, "[traceback: %.2lfms] ", time);
    time= pop_time();
    fprintf(stdout, "done. [%.2lfms] [cost: %d] ", time, min_cost);
    print_pop_memtrack();
    fprintf(stdout, " ");
    print_current_mem();
    fprintf(stdout, "\n");
    fflush(stdout);
    // list a solution
    list_solution(g);
    graph_free(g);
#endif

    FREE(d);
    FREE(f_v); 
    FREE(visit);
    FREE(p);
    FREE(b_v);
#ifdef TRACK_BANDWIDTH
    FREE(heap_ops);
#endif

    return min_cost;
}

/******************************************************* Program entry point. */

int main(int argc, char **argv)
{
#ifdef TRACK_RESOURCES
    push_time();
	push_memtrack();

    fprintf(stdout, "invoked as:");
    for(index_t f = 0; f < argc; f++) 
        fprintf(stdout, " %s", argv[f]);
    fprintf(stdout, "\n");
#endif

    if(!strcmp(argv[1], "-h")) {
        fprintf(stdout, "Usage: %s -s <seed> <in-file>\n\n", argv[0]);
        return 0;
    }

    FILE *in = NULL;
    if(argc < 4) {
        //ERROR("Insufficient arguments");
        // hack: read graph from standard input
        in = stdin;
    }
    else {
        char *filename = argv[3];
        in = fopen(filename, "r");
        if(in == NULL)
            ERROR("unable to open file '%s'", filename);
    }

    // read input graph
    graph_t *g = graph_load(in);
    // build root query
    steinerq_t *root = root_build(g);
    // release graph memory
    graph_free(g);
    // execute the algorithm
    erickson_monma_veinott(root);
    // release query memory
    steinerq_free(root);

#ifdef TRACK_RESOURCES
    double time = pop_time();
    fprintf(stdout, "grand total [%.2lfms] ", time);
    print_pop_memtrack();
    fprintf(stdout, "\n");
    fprintf(stdout, "host: %s\n", sysdep_hostname());
    fprintf(stdout, "build: %s, %s, %s\n",
                    "edge-linear kernel",
                    "single thread",
                    "binary heap"
            );
    fprintf(stdout,
            "compiler: gcc %d.%d.%d\n",
            (index_t)__GNUC__,
            (index_t)__GNUC_MINOR__,
            (index_t)__GNUC_PATCHLEVEL__);
    fflush(stdout);
    assert(malloc_balance == 0);
    assert(memtrack_stack_top < 0);
#endif
    return 0;
}
