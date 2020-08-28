import numpy as np
import msprime
import sys


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

    assert(last_tmrca != current_tmrca)

    update_matrices_and_output_ibd_for_sample_pair.longest = max(
        update_matrices_and_output_ibd_for_sample_pair.longest,
        IBD_seg_length)

    if IBD_seg_length > ibd_threshold_in_bp:
        print("{id1}\t{hap1}\t{id2}\t{hap2}\t{chrom}\t"
              "{start:0.2f}\t{end:0.2f}\t{cM:0.2f}\t{tmrca:0.2f}".format(
                  id1=row, hap1=1, id2=col, hap2=1,
                  chrom=chrom, start=last_break_pos,
                  end=interval_start_pos,
                  cM=IBD_seg_length / bp_per_cM,
                  tmrca=last_tmrca))
        update_matrices_and_output_ibd_for_sample_pair.num_ibd += 1


def output_last_ibds_for_all_pairs(
        tmrca_matrix, breakpos_matrix, num_samples, tree_end,
        chrom, bp_per_cM, ibd_threshold_in_bp):

    for col in range(0, num_samples - 1):
        for row in range(col + 1, num_samples):
            tmrca = tmrca_matrix[row, col]
            last_break_pos = breakpos_matrix[row, col]
            IBD_seg_length = tree_end - last_break_pos
            if IBD_seg_length > ibd_threshold_in_bp:
                print("{id1}\t{hap1}\t{id2}\t{hap2}\t{chrom}\t"
                      "{start:0.2f}\t{end:0.2f}\t{cM:0.2f}\t"
                      "{tmrca:0.2f}".format(
                          id1=row, hap1=1, id2=col, hap2=1,
                          chrom=chrom, start=last_break_pos,
                          end=tree_end,
                          cM=IBD_seg_length / bp_per_cM,
                          tmrca=tmrca))


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
        tree, M, RN, tmrca_matrix, breakpos_matrix,
        chrom, bp_per_cM, ibd_threshold_in_bp):

    interval_start_pos = tree.interval[0]
    new_leaves = find_tree_leaves_of_any_node(tree, M)
    run_node = M
    num_pairs_updated = 0
    while run_node != RN:
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
    # for edges_in
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
        M = N
        RN = N
        pat = "Pattern 1: two-edge, moving root itself"

    # condition 2: three-edge patterns
    elif num_edges == 3:
        two_children_Nodes = [parent for parent in in_pc_map
                              if len(in_pc_map[parent]) == 2]
        two_children_Nodes_out = [parent for parent in out_pc_map
                                  if len(out_pc_map[parent]) == 2]
        # find N node:
        Ns = [node for node in two_children_Nodes
              if node in in_cp_map or node == tree.root]
        assert(len(Ns) == 1)
        N = Ns[0]
        Pn = tree.parent(N)

        # this is specific for pattern 2B
        if len(two_children_Nodes_out) == 1:
            R = two_children_Nodes_out[0]
        else:
            R = None

        # # condition 2A: N node is on the edge of the prev tree root;
        # thus prev tree root is larger number of N node's children of the new
        # tree, the M node the smaller number

        if N == tree.root:
            M = min(in_pc_map[N])
            RN = N
            pat = "Pattern 2A: three-edge, moving to the root"

        # # condition 2B: N node is along the edge of R node in old tree; N
        # node in new tree  and R in old tree have same children

        elif N != tree.root and R is not None and R in out_cp_map \
                and out_cp_map[R] == in_cp_map[N]:
            assert(len(in_pc_map[N]) == 2)
            M = in_pc_map[N][1]
            RN = N  # no need to update for nodes older than N, thus RN=N
            pat = "Pattern 2B: three-edge, moving along an edge"

        # # condition 2C: similar to four-edge pattern (only diff: moving node
        # is a child node of the root
        else:
            B_list = list(set(out_pc_map[Pn]).intersection(set(in_pc_map[N])))
            assert(len(B_list) == 1)
            B = B_list[0]
            # M node is N's other child, B's sibling
            N_children = in_pc_map[N]
            assert(len(N_children) == 2)
            if B == N_children[0]:
                M = N_children[1]
            else:
                M = N_children[0]
                assert(B == N_children[1])
            # R is the root of old tree, movement of M destroy one child;
            # R's other child and M's sibling of the old treeis the root of new
            # tree, RN is root of the new tree
            RN = tree.root
            pat = "Pattern 2C: three-edge, moving a root child to an edge"

    # condition 3: four-edge pattern
    else:  #
        # N node
        two_children_Nodes = [parent for parent in in_pc_map
                              if len(in_pc_map[parent]) == 2]
        N_list = [node for node in two_children_Nodes if node in in_cp_map]
        assert(len(N_list) == 1)
        N = N_list[0]
        # B node
        Pn = tree.parent(N)
        B_list = list(set(out_pc_map[Pn]).intersection(set(in_pc_map[N])))
        assert(len(B_list) == 1)
        B = B_list[0]
        # M node
        N_children = in_pc_map[N]
        assert(len(N_children) == 2)
        if B == N_children[0]:
            M = N_children[1]
        else:
            M = N_children[0]
            assert(B == N_children[1])
        # RN node
        R = out_cp_map[M]  # R is parent of M in old tree
        Rp = out_cp_map[R]
        RN = tree.mrca(Rp, N)
        pat = "Pattern 3: four-edge, general"

    # print("RN: {}. Is root? {}".format(RN, RN == tree.root))

    return M, N, RN, pat


def test_tsibd_incremental_vs_complete():

    # some basic options
    chrom = 1
    bp_per_cM = 1000000
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

    # initialize function attributes
    update_matrices_and_output_ibd_for_sample_pair.longest = 0
    update_matrices_and_output_ibd_for_sample_pair.num_ibd = 0

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
        # if int(interval[0]) / 10000 == int(interval[1]) / 10000:
        #     continue

        print("tree: %d, position: %f, longest: %f, num_ibd: %f" % (
            i, interval[1],
            update_matrices_and_output_ibd_for_sample_pair.longest,
            update_matrices_and_output_ibd_for_sample_pair.num_ibd
        ))

        M, N, RN, pat = find_edge_diff_pattern(diff, tree)

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
                    tree, M, RN, tmrca_matrix, breakpos_matrix,
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

    output_last_ibds_for_all_pairs(
        tmrca_matrix, breakpos_matrix, num_samples,
        interval[1],  # last interval
        chrom, bp_per_cM, ibd_threshold_in_bp)

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


def test_tsibd():

    # some basic options
    chrom = 1
    bp_per_cM = 1000000
    ibd_threshold_in_bp = 2000000

    # simulation
    tree_sequence = msprime.simulate(
        sample_size=100, Ne=10000, length=30000000, recombination_rate=1e-8)

    print("Done with simulation!", file=sys.stderr)

    # make matrics
    num_samples = tree_sequence.num_samples
    tmrca_matrix = np.zeros(shape=(num_samples, num_samples))
    breakpos_matrix = np.zeros(shape=(num_samples, num_samples))

    # initialize function attributes
    update_matrices_and_output_ibd_for_sample_pair.longest = 0
    update_matrices_and_output_ibd_for_sample_pair.num_ibd = 0

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

    interval = None
    for i in range(1, tree_sequence.num_trees):

        # get the edge_diff and tree
        interval = tree.interval
        diff = next(diffs)
        tree = next(trees)

        if int(interval[0]) / 10000 == int(interval[1]) / 10000:
            continue

        print("tree: %d, position: %f, longest: %f, num_ibd: %f" % (
            i, interval[1],
            update_matrices_and_output_ibd_for_sample_pair.longest,
            update_matrices_and_output_ibd_for_sample_pair.num_ibd
        ), file=sys.stderr)

        M, N, RN, pat = find_edge_diff_pattern(diff, tree)

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
                    tree, M, RN, tmrca_matrix, breakpos_matrix,
                    chrom, bp_per_cM, ibd_threshold_in_bp)

    output_last_ibds_for_all_pairs(
        tmrca_matrix, breakpos_matrix, num_samples,
        interval[1],  # last interval
        chrom, bp_per_cM, ibd_threshold_in_bp)

    print("", file=sys.stderr)
    print("-----------------", file=sys.stderr)
    print("sample_size: %d" % num_samples, file=sys.stderr)
    print("num_pairs_updated_incremental: %d" % num_pairs_updated_incremental,
          file=sys.stderr)


if __name__ == "__main__":

    """
    comparing incremental vs complete
    """

    # ratio_list = []
    # num_test = 1
    # for i in range(num_test):
    #     ratio = test_tsibd_incremental_vs_complete()

    #     if ratio is not None:
    #         ratio_list.append(ratio)

    # print("\n ============ ")
    # print("num_test", num_test)
    # print("num valid ratios", len(ratio_list))
    # print("Ratios: ")
    # print("     mean: %g" % np.mean(np.array(ratio_list)))
    # print("     std: %g" % np.std(np.array(ratio_list)))

    """
    just run incremental
    """
    test_tsibd()
