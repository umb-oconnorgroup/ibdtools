import numpy as np
import sys
import msprime


def find_tree_leaves_of_any_node(tree, u):
    if tree.is_sample(u):
        leaves = [u]
    else:
        leaves = list(tree.leaves(u))
    return leaves


def update_matrices_and_output_ibd_for_sample_pair(
        tmrca_matrix, breakpos_matrix, id1, id2, current_tmrca,
        interval_start_pos, chrom, bp_per_cM,
        ibd_threshold_in_bp):

    if id1 < id2:
        row, col = id2, id1
    else:
        row, col = id1, id2

    """
    TODO: codes deal with special case like the last tree. last tree all should
    breaks at the end of the interval, but current implementation igored them..
    """

    last_tmrca = tmrca_matrix[row, col]
    tmrca_matrix[row, col] = int(current_tmrca)
    last_break_pos = breakpos_matrix[row, col]
    breakpos_matrix[row, col] = interval_start_pos
    IBD_seg_length = interval_start_pos - last_break_pos

    update_matrices_and_output_ibd_for_sample_pair.longest = max(
        update_matrices_and_output_ibd_for_sample_pair.longest,
        IBD_seg_length)

    if IBD_seg_length > ibd_threshold_in_bp:
        print("{id1}\t{hap1}\t{id2}\t{hap2}\t{chr}\t"
              "{start}\t{end}]\t{cM}\t{tmrca}".format(
                  id1=row, hap1=1, id2=col, hap2=1,
                  chr=1, start=last_break_pos,
                  end=interval_start_pos,
                  cM=IBD_seg_length / bp_per_cM,
                  tmrca=last_tmrca))
        update_matrices_and_output_ibd_for_sample_pair.num_ibd += 1

update_matrices_and_output_ibd_for_sample_pair.longest = 0
update_matrices_and_output_ibd_for_sample_pair.num_ibd = 0


def update_matrices_due_to_new_ancestral_node(
        tree, u, tmrca_matrix, breakpos_matrix,
        chrom, bp_per_cM, ibd_threshold_in_bp):
    assert(not tree.is_sample(u))
    current_tmrca = tree.time(u)
    interval_start_pos = tree.interval[0]
    left_child = tree.left_child(u)
    right_child = tree.right_child(u)
    left_child_leaves = find_tree_leaves_of_any_node(tree, left_child)
    right_child_leaves = find_tree_leaves_of_any_node(tree, right_child)
    for id1 in left_child_leaves:
        for id2 in right_child_leaves:
            update_matrices_and_output_ibd_for_sample_pair(
                tmrca_matrix, breakpos_matrix, id1, id2,
                current_tmrca, interval_start_pos, chrom,
                bp_per_cM, ibd_threshold_in_bp)
    num_pairs_updated = len(left_child_leaves) * len(right_child_leaves)
    return num_pairs_updated


def update_matries_due_to_new_leaves(
        tree, src_node, tmrca_matrix, breakpos_matrix,
        chrom, bp_per_cM, ibd_threshold_in_bp):

    interval_start_pos = tree.interval[0]
    new_leaves = find_tree_leaves_of_any_node(tree, src_node)
    run_node = src_node
    num_pairs_updated = 0
    while run_node != tree.root:
        parent = tree.parent(run_node)
        current_tmrca = tree.time(parent)
        if run_node == tree.left_child(parent):
            sib = tree.right_child(parent)
        else:
            sib = tree.left_child(parent)
        other_side_leaves = find_tree_leaves_of_any_node(tree, sib)
        for id1 in new_leaves:
            for id2 in other_side_leaves:
                update_matrices_and_output_ibd_for_sample_pair(
                    tmrca_matrix, breakpos_matrix, id1, id2,
                    current_tmrca, interval_start_pos, chrom,
                    bp_per_cM, ibd_threshold_in_bp)
        num_pairs_updated += len(new_leaves) * len(other_side_leaves)
        run_node = parent

    return num_pairs_updated


def update_tmrca_matrix_due_to_new_tree(
        tree, tmrca_matrix, breakpos_matrix,
        chrom, bp_per_cM, ibd_threshold_in_bp):

    num_pairs_updated = 0
    for u in tree.nodes():
        # access internal nodes
        if not tree.is_sample(u):
            num_pairs_updated += update_matrices_due_to_new_ancestral_node(
                    tree, u, tmrca_matrix, breakpos_matrix,
                    chrom, bp_per_cM, ibd_threshold_in_bp)
    return num_pairs_updated


def print_diff(edge_diff, j):
    edges_out = edge_diff[1]
    edges_in = edge_diff[2]

    print("Edge_diffs: between tree %d and tree %d" % (j - 1, j))
    print("    edges out", end=": ")
    for edge in edges_out:
        print("%d -> %d;" % (edge.child, edge.parent), end=' ')
    print("")
    print("    edges in", end=": ")
    for edge in edges_in:
        print("%d -> %d;" % (edge.child, edge.parent), end=' ')
    print("")


def print_tree(tree, j):
    print("                                    Tree: %d" % j)
    print(tree.draw(format="unicode"))


def find_edge_diff_pattern(diff, tree):
    """
    Most complicated algorithm here: is to find the key nodes from the
    edge_diffs record
    """

    # build parent to child and child to parent mapping for edges_out and
    # edges_int
    edges_out = diff[1]
    edges_in = diff[2]
    in_cp_map = {}
    in_pc_map = {}
    out_cp_map = {}
    out_pc_map = {}

    for edge in edges_in:
        in_cp_map[edge.child] = edge.parent
        if edge.parent not in in_pc_map:
            in_pc_map[edge.parent] = [edge.child]
        else:
            in_pc_map[edge.parent].append(edge.child)

    for edge in edges_out:
        out_cp_map[edge.child] = edge.parent
        if edge.parent not in out_pc_map:
            out_pc_map[edge.parent] = [edge.child]
        else:
            out_pc_map[edge.parent].append(edge.child)

    num_edges = len(diff[1])

    # ========== Find the N and M node ===============
    # condition 1: two-edge pattern: movement tree root
    if num_edges == 2:  # only two edges
        N = edges_in[0].parent
        M = edges_in[0].parent
        pat = "Pattern 1: two-edge, moving root itself"

    # condition 2: three-edge patterns
    elif num_edges == 3:
        # first find N node:
        two_children_Nodes = [parent for parent in in_pc_map
                              if len(in_pc_map[parent]) == 2]
        two_children_Nodes_out = [parent for parent in out_pc_map
                                  if len(out_pc_map[parent]) == 2]
        if len(two_children_Nodes_out) == 1:
            R_2B = two_children_Nodes_out[0]
        else:
            R_2B = None

        Ns = [node for node in two_children_Nodes
              if node in in_cp_map or node == tree.root]
        assert(len(Ns) == 1)
        N = Ns[0]
        Pn = tree.parent(N)

        # # condition 2A: N node is on the edge of the prev tree root
        # which is larger number of N node's children of the new tree, the M
        # node the smaller number
        if N == tree.root:
            M = min(in_pc_map[N])
            pat = "Pattern 2A: three-edge, moving to the root"

        # # condition 2B: moving node just move along the original edge
        # R node in old tree and N node in new tree have same children
        elif N != tree.root and R_2B is not None and R_2B in out_cp_map \
                and out_cp_map[R_2B] == in_cp_map[N]:
            assert(len(in_pc_map[N]) == 2)
            M = in_pc_map[N][1]
            pat = "Pattern 2B: three-edge, moving along an edge"

        # # condition 2C: similar to four-edge pattern (only diff: moving node
        # is a child node of the root
        else:
            B_proposal = set(out_pc_map[Pn]).intersection(set(in_pc_map[N]))
            B_proposal = list(B_proposal)
            if(len(B_proposal) != 1):
                print("something is wrong with B")
                B = -1
            else:
                B = B_proposal[0]
            # M node
            N_children = in_pc_map[N]
            assert(len(N_children) == 2)
            if B == N_children[0]:
                M = N_children[1]
            else:
                M = N_children[0]
                # assert(B == N_children[1])
                if B != N_children[1]:
                    print("some things is wrong with M")
            pat = "Pattern 2C: three-edge, moving a root child to an edge"

    # condition 3: four-edge pattern
    else:  #
        # N node
        two_children_Nodes = [parent for parent in in_pc_map
                              if len(in_pc_map[parent]) == 2]
        Ns = [node for node in two_children_Nodes
              if node in in_cp_map]
        assert(len(Ns) == 1)
        N = Ns[0]
        # B node
        Pn = tree.parent(N)
        B_proposal = set(out_pc_map[Pn]).intersection(set(in_pc_map[N]))
        B_proposal = list(B_proposal)
        B = B_proposal[0]
        if(len(B_proposal) != 1):
            print("something is wrong with B")
        # M node
        N_children = in_pc_map[N]
        assert(len(N_children) == 2)
        if B == N_children[0]:
            M = N_children[1]
        else:
            M = N_children[0]
            # assert(B == N_children[1])
            if B != N_children[1]:
                print("some things is wrong with M")
        pat = "Pattern 3: four-edge, general"

    return M, N, pat


def test_tmrca_update_incremental_vs_complete():

    # some basic options
    chrom = 1
    bp_per_cM = 1000000,
    ibd_threshold_in_bp = 2000000

    # simulation
    tree_sequence = msprime.simulate(
        sample_size=100, Ne=10000, length=30000000, recombination_rate=1e-8)

    print("Done with simulation!")

    # make matrics
    num_samples = tree_sequence.num_samples
    tmrca_matrix = np.zeros(shape=(num_samples, num_samples))
    breakpos_matrix = np.zeros(shape=(num_samples, num_samples))
    breakpos_matrix_complete = np.zeros(shape=(num_samples, num_samples))

    # initialize the tmrca matrix from first tree
    first_tree = tree_sequence.first()
    update_tmrca_matrix_due_to_new_tree(
        first_tree, tmrca_matrix, breakpos_matrix,
        chrom, bp_per_cM, ibd_threshold_in_bp)
    # print(tmrca_matrix)

    # updates matrix for following trees using tree_diffs
    trees = tree_sequence.trees()
    diffs = tree_sequence.edge_diffs()

    # skip the first tree in the for loop
    tree = next(trees)
    diff = next(diffs)

    # num_pair_updated_coutners initialization
    num_pairs_updated_incremental = 0
    num_pairs_updated_complete = 0

    matrix_all_equal_list = []

    # print_tree(tree, 0)
    for i in range(1, tree_sequence.num_trees):

        # get the edge_diff and tree
        interval = tree.interval
        diff = next(diffs)
        tree = next(trees)

        # This is to skip tree every 10k: the disadvantage is that this makes
        # the complete update and incremental update different.
        if int(interval[0]) / 10000 == int(interval[1]) / 10000:
            continue

        print("tree: %d, position: %f, longest: %f, num_ibd: %f" % (
            i, interval[1],
            update_matrices_and_output_ibd_for_sample_pair.longest,
            update_matrices_and_output_ibd_for_sample_pair.num_ibd
        ))

        M, N, pat = find_edge_diff_pattern(diff, tree)

        # NOTE: Here are the key differences

        # 1. update exisiting trmca_matrix (incremental)
        if M == tree.root:
            num_pairs_updated_incremental += \
                update_matrices_due_to_new_ancestral_node(
                    tree, M, tmrca_matrix, breakpos_matrix,
                    chrom, bp_per_cM, ibd_threshold_in_bp)
        else:
            num_pairs_updated_incremental += \
                update_matries_due_to_new_leaves(
                    tree, M, tmrca_matrix, breakpos_matrix,
                    chrom, bp_per_cM, ibd_threshold_in_bp)

        # 2. make an brand new matrix from the current tree (complete)
        tmrca_matrix_temp = np.zeros(shape=(num_samples, num_samples))
        num_pairs_updated_complete += \
            update_tmrca_matrix_due_to_new_tree(
                tree, tmrca_matrix_temp, breakpos_matrix_complete,
                chrom, bp_per_cM, ibd_threshold_in_bp)

        # compare two matrices

        # NOTE: skipping trees will make the complete update and incremental
        # update different.
        matrix_all_equal_list.append(np.all(tmrca_matrix == tmrca_matrix_temp))
        all_equal = np.all(np.array(matrix_all_equal_list))

        """ FOR DEBUGGING
        if not all_equal:
            print("matrix not equal")

            tree.prev()
            print_tree(tree, i-1)
            tree.next()
            print("               |")
            print("               |")
            print("               v")
            print_diff(diff, i)
            print("")
            print(pat)
            print(tmrca_matrix)
            print(tmrca_matrix_temp)
            print("Is incremental TMRCA matrix == complete TMRCA matrix?")
            print(tmrca_matrix == tmrca_matrix_temp)
            print("** N Node: %d, M Node: %d **" % (N, M))
            print("               |")
            print("               |")
            print("               v")
            print_tree(tree, i)
         """

    print("")
    print("-----------------")
    print("sample_size", num_samples)
    print("two methods generate the same matrix for each marginal tree?",
          all_equal)
    print("num_pairs_updated_complete", num_pairs_updated_complete)
    print("num_pairs_updated_incremental", num_pairs_updated_incremental)
    if num_pairs_updated_complete != 0:
        ratio = num_pairs_updated_incremental / num_pairs_updated_complete
        print("incremental to complete ratio: ", ratio)
    else:
        ratio = None
    return ratio


if __name__ == "__main__":
    ratio_list = []
    num_test = 1
    for i in range(num_test):
        ratio = test_tmrca_update_incremental_vs_complete()
        if ratio is not None:
            ratio_list.append(ratio)

    print("\n ============ ")
    print("num_test", num_test)
    print("num valid ratios", len(ratio_list))
    print("Ratios: ")
    print("     mean: %g" % np.mean(np.array(ratio_list)))
    print("     std: %g" % np.std(np.array(ratio_list)))
