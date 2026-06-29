#include <stdio.h>
#include <stdlib.h>
#include <sys/queue.h>
#include <stdbool.h>
#include <dotgeno.h>

typedef struct {
	snp_data* elems;
	size_t length;
} snp_objs;

typedef struct {
	struct idx_head* ref_idx;
	struct idx_head* elm_idx;  // come up with better name
} idx_pair;

typedef struct {
	idx_pair* ip;
	size_t length;
} idx_intersect;

snp_objs init_snp_objs(size_t length) {
	snp_objs snps;
	snps.length = length;
	snps.elems = (snp_data*)malloc(sizeof(snp_data) * snps.length);
	return snps;
}

idx_intersect init_idx_intersect(snp_objs* snps) {
	idx_intersect isct;
	isct.length = snps->length - 1;
	isct.ip = (idx_pair*)malloc(sizeof(idx_pair) * isct.length);
	return isct;
}

idx_pair get_idx_pair(snp_data* ref, snp_data* cmp) {
	idx_pair pair;
	pair.ref_idx = malloc(sizeof *pair.ref_idx);
	pair.elm_idx = malloc(sizeof *pair.elm_idx);

	STAILQ_INIT(pair.ref_idx);
	STAILQ_INIT(pair.elm_idx);
	intersect_snp_data(ref, cmp, pair.ref_idx, pair.elm_idx);
	return pair;
}

size_t get_max_index(struct idx_node** arr, size_t length) {
	size_t max = 0;
	size_t max_idx = 0;
	for(size_t i = 0; i < length; i++) {
		if(arr[i]->idx >= max) {
			max = arr[i]->idx;
			max_idx = i;
		}
	}
	return max_idx;
}

bool all_equal(struct idx_node** arr, size_t length) {
	size_t first = arr[0]->idx;
	for(size_t i = 1; i < length; i++) {
		if(first != arr[i]->idx) {
			return false;
		}
	}
	return true; 
}


void print_isct(idx_intersect* isct) {
	for(size_t i = 0; i < isct->length; i++) {
		printf("Index %u:\n", i);
		printf("\tIDX:");
		struct idx_node* in;
		STAILQ_FOREACH(in, isct->ip[i].ref_idx, nodes) {
			printf(" %u", in->idx);
		}
		printf("\n");
	}
}

void make_intersection(idx_intersect* isct) {
	struct idx_node** cur_elems = (struct idx_node**)malloc(isct->length * sizeof(struct idx_node*));
	struct idx_node** cur_elems_act = (struct idx_node**)malloc(isct->length * sizeof(struct idx_node*));
	for(size_t i = 0; i < isct->length; i++) {
		cur_elems[i] = STAILQ_FIRST(isct->ip[i].ref_idx);
		cur_elems_act[i] = STAILQ_FIRST(isct->ip[i].elm_idx);
	}
	print_isct(isct);
	while(1) {
		bool end_while = false;
		if(all_equal(cur_elems, isct->length)) {
			for(size_t i = 0; i < isct->length; i++) {
				struct idx_node* n_a = STAILQ_NEXT(cur_elems[i], nodes);
				struct idx_node* n_b = STAILQ_NEXT(cur_elems_act[i], nodes);
				if(n_a == NULL) {
					for(size_t j = i + 1; j < isct->length; j++) {
						if(j == i) { continue; }
						struct idx_node* nxt_a = STAILQ_NEXT(cur_elems[j], nodes);
						struct idx_node* nxt_b = STAILQ_NEXT(cur_elems_act[j], nodes);
						while(nxt_a) {
							STAILQ_REMOVE(isct->ip[j].ref_idx, nxt_a, idx_node, nodes);
							STAILQ_REMOVE(isct->ip[j].elm_idx, nxt_b, idx_node, nodes);
							nxt_a = STAILQ_NEXT(cur_elems[j], nodes);
							nxt_b = STAILQ_NEXT(cur_elems_act[j], nodes);
						}
					}
					for(size_t j = 0; j < i; j++) {
						if(j == i) { continue; }
						struct idx_node* nxt_a = cur_elems[j];
						struct idx_node* nxt_b = cur_elems_act[j];
						while(nxt_a != NULL) {
							STAILQ_REMOVE(isct->ip[j].ref_idx, nxt_a, idx_node, nodes);
							STAILQ_REMOVE(isct->ip[j].elm_idx, nxt_b, idx_node, nodes);
							nxt_a = STAILQ_NEXT(cur_elems[j], nodes);
							nxt_b = STAILQ_NEXT(cur_elems_act[j], nodes);
						}
					}
					end_while = true;
					break;
				}
				cur_elems[i] = n_a;
				cur_elems_act[i] = n_b;
			}
		} else {
			size_t max_i = get_max_index(cur_elems, isct->length);
			for(size_t i = 0; i < isct->length; i++) {
				if(i == max_i) { continue; }
				// remove and move to the next one
				struct idx_node* n_a = STAILQ_NEXT(cur_elems[i], nodes);
				struct idx_node* n_b = STAILQ_NEXT(cur_elems_act[i], nodes);
				STAILQ_REMOVE(isct->ip[i].ref_idx, cur_elems[i], idx_node, nodes);
				STAILQ_REMOVE(isct->ip[i].elm_idx, cur_elems_act[i], idx_node, nodes);
				free(cur_elems[i]);
				free(cur_elems_act[i]);
				if(n_a == NULL) {
					for(int j = 0; j < isct->length; j++) {
						free_idx_list(isct->ip[j].ref_idx);
						free_idx_list(isct->ip[j].elm_idx);
					}
					end_while = true;
					break;
				}
				cur_elems[i] = n_a;
				cur_elems_act[i] = n_b;
			}
		}
		print_isct(isct);
		if(end_while) {
			free(cur_elems);
			free(cur_elems_act);
			break;
		}
	}
}

int main(int argc, char* argv[]) {
	snp_objs snps = init_snp_objs(argc - 1);
	for(size_t i = 0; i < snps.length; i++) {
		snps.elems[i] = read_snp_file(argv[i+1]);
	}
	idx_intersect isct = init_idx_intersect(&snps);
	for(size_t i = 0; i < isct.length; i++) {
		isct.ip[i] = get_idx_pair(&snps.elems[0], &snps.elems[i+1]);
	}
	make_intersection(&isct);
}
