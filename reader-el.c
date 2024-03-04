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

typedef int index_t;

/********************************************************** Global constants. */

#define MAX_DISTANCE ((index_t)0x7FFFFFFF)
#define MATH_INF ((index_t)0x7FFFFFFF)
#define UNDEFINED -1

/************************************************************* Common macros. */

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

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

    return g;
}

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

    for(index_t u = 0; u < n; u++)
        pos[u] = 0;

    index_t *e = g->edges;
    for(index_t j = 0; j < 3*m; j+=3)
    {
        pos[e[j]]+=2;
        pos[e[j+1]]+=2;
    }

    pos[n] = (2*n);
    index_t run = prefixsum(n+1, pos, 1);
    assert(run == ((n+1)+(4*m)+(2*n)));

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

    for(index_t u = 0; u < k; u++)
        kk[u] = tt[u];

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

} bheap_t;

bheap_t * bh_alloc(index_t n)
{
    bheap_t *h = (bheap_t *) malloc(sizeof(bheap_t));
    h->max_n = n;
    h->n = 0; 
    h->a = (bheap_item_t *) malloc((n+1)*sizeof(bheap_item_t));
    h->p = (index_t *) malloc(n*sizeof(index_t));


    return h;
}

void bh_free(bheap_t *h)
{
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

    while(k <= q)
    {
        bheap_item_t z = h->a[k];
        if(k < q)
        {
            if(z.key > h->a[k + 1].key) z = h->a[++k];
        }
        if(y.key <= z.key) break;
        h->a[j] = z;
        h->p[z.item] = j;
        j = k;
        k = 2 * j;
    }

    h->a[j] = y;
    h->p[y.item] = j;
}

bheap_item_t bh_min(bheap_t *h)
{
    return (bheap_item_t) h->a[1];
}

static void bh_insert(bheap_t *h, index_t item, index_t key)
{
    index_t i = ++(h->n);
    while(i >= 2)
    {
        index_t j = i / 2;
        bheap_item_t y = h->a[j];

        if(key >= y.key) break;

        h->a[i] = y;
        h->p[y.item] = i;
        i = j;
    }

    h->a[i].item = item;
    h->a[i].key = key;
    h->p[item] = i;
}

static void bh_delete(bheap_t *h, index_t item)
{
    index_t n = --(h->n);
    index_t p = h->p[item];

    if(p <= n)
    {
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
}

static void bh_decrease_key(bheap_t *h, index_t item, index_t new_key)
{
    index_t i = h->p[item];
    while(i >= 2)
    {
        index_t j = i / 2;
        bheap_item_t y = h->a[j];
        if(new_key >= y.key) break;

        h->a[i] = y;
        h->p[y.item] = i;
        i = j;
    }

    h->a[i].item = item;
    h->a[i].key = new_key;
    h->p[item] = i;
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

void backtrack(index_t n, index_t k, index_t v, index_t X, index_t *kk, index_t *b_v)
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
        backtrack(n, k, u, Xd, kk, b_v);
    }
    else
    {
        index_t Xd   = b_v[i_X+1];
        index_t X_Xd = (X & ~Xd);
        if(X == Xd)
            return;
        backtrack(n, k, u, Xd, kk, b_v);
        backtrack(n, k, u, X_Xd, kk, b_v);
    }
}

void build_tree(index_t n, index_t k, index_t cost, index_t *kk, index_t *b_v)
{
    index_t c = k-1;
    index_t C = (1<<c)-1;
    index_t q = kk[k-1];

#ifdef TRACK_RESOURCES
    g->n = n;
#else
    fprintf(stdout, "VALUE %d\n", cost);
#endif

    backtrack(n, k, q, C, kk, b_v);
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
                   index_t *b_v)
{
    // initialisation
    index_t c = k-1;
    index_t q = kk[k-1];
    index_t C = (1<<c)-1;

    for(index_t t = 0; t < k; t++) 
    {    
        dijkstra(n+1, m, pos, adj, kk[t], d, visit, p);
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
        // generating all subsets of size `m`
        for(index_t X = (1<<m)-1;
            X < (1<<c);
            z = X|(X-1), X = (z+1)|(((~z & -~z)-1) >> (__builtin_ctz(X) + 1))) // cCm
        {
            index_t *f_X  = f_v + FV_INDEX(0, n, k, X);
            index_t *b_X  = b_v + BV_INDEX(0, n, k, X);
            // generating proper subsets of X 
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

            dijkstra(n+1, m+n, pos, adj, s, d, visit, p);
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

    dijkstra(n+1, m+n, pos, adj, s, d, visit, p);
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

        dijkstra(n, m, root->pos, root->adj, u, d, visit, p);

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

    // initialisation
    for(index_t i = 0; i < (index_t)(n*(1<<c)); i++)
        f_v[i] = MATH_INF;

    for(index_t i = 0; i < (index_t)(2*n*(1<<c)); i+=2) 
    {
        b_v[i]   = -1;
        b_v[i+1] = 0;
    }

    // call kernel
    min_cost = emv_kernel(n, m, k, root->pos, root->adj, kk, d, p, visit, f_v, b_v);
    // build a Steiner tree
    build_tree(n, k, min_cost, kk, b_v);

    FREE(d);
    FREE(f_v); 
    FREE(visit);
    FREE(p);
    FREE(b_v);

    return min_cost;
}

int main(int argc, char **argv)
{
    if(!strcmp(argv[1], "-h")) {
        fprintf(stdout, "Usage: %s -s <seed> <in-file>\n\n", argv[0]);
        return 0;
    }

    FILE *in = NULL;
    if(argc < 4) {
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

    return 0;
}
