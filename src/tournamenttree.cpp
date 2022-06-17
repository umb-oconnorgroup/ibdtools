#include "common.hpp"

template <typename T>
TournamentTree<T>::TournamentTree(int k_, T max_val) : max_val(max_val), k(k_)
{
    exit_on_false(
        k >= 1, "the number of node needs to be greater than 1", __FILE__, __LINE__);

    for (N = 1; (1UL << N) < (size_t) k; N++)
        ;
    // std::cout << "k = " << k << " N: " << N << " 2^n = " << (1 << N) << '\n';
    size_t sz = 0;
    for (size_t layer = 0; layer <= N; layer++)
        sz += (1 << layer);
    // std::cout << "size = " << sz << '\n';
    node_vec.resize(sz);

    size_max = std::numeric_limits<size_t>::max();
}

/*
 *                       0
 *           1                       2
 *     3           4           5           6
 *  7     8     9     10    11    12    13    14
 */
// @ winner_id receives the next winner. (output)
// @ initial_val: is a vector values to initial the K slots.
// @ return the the winner's value.
template <typename T>
[[nodiscard]] T &
TournamentTree<T>::init_run(std::vector<T> initial_vals, size_t &winner_id)
{
    exit_on_false(k == initial_vals.size(), "", __FILE__, __LINE__);
    std::vector<size_t> layer_start_vec;
    for (size_t layer = 0, count = 0; layer <= N; count += (1 << layer), layer++) {
        layer_start_vec.push_back(count);
        // std::cout << "layer: " << layer << "  start_id: " << count << '\n';
    }
    layer_N_start = layer_start_vec.back();

    for (size_t layer = N; layer != size_max; layer--) {
        // std::cout << "for layer: " << layer << '\n';
        for (size_t i = 0; i < (1UL << layer); i++) {
            size_t this_node_id = layer_start_vec[layer] + i;
            Node &node = node_vec[this_node_id];
            // initialize parent
            if (layer > 0) {
                size_t prev_layer_start = layer_start_vec[layer - 1];
                node.parent = prev_layer_start + i / 2;
            } else {
                node.parent = size_max;
            }
            // initialize leaf and initialize value_id
            if (layer == N) {
                node.leaf = i;
                if (i < initial_vals.size())
                    node.val = initial_vals[i];
                else {
                    node.val = max_val;
                }
            } else {
                size_t next_layer_start = layer_start_vec[layer + 1];
                size_t left_child_node_id = next_layer_start + 2 * i;
                size_t right_child_node_id = left_child_node_id + 1;
                Node &left = node_vec[left_child_node_id];
                Node &right = node_vec[right_child_node_id];
                if (left.val < right.val) {
                    node.leaf = left.leaf;
                    node.val = left.val;
                } else {
                    node.leaf = right.leaf;
                    node.val = right.val;
                }
            }

            // std::cout << "node: " << this_node_id << "\t i: " << i
            //           << "\t leaf: " << node.leaf << "\t val: " << node.val
            //           << "\t parent: " << node.parent << '\n';
        }
    }
    winner_id = root().leaf;
    return root().val;
}

// @ winner_id receives the next winner. (output)
// @ T val is the input number to replace the winner's slot.
// @ return the the winner's value.
template <typename T>
T &
TournamentTree<T>::replace_run(T val, size_t &winner_id)
{
    size_t this_node_id = layer_N_start + root().leaf;
    // size_t parent_id, sister_id;

    // update this node's value
    node_vec[this_node_id].val = val;

    // update parents
    for (size_t layer = N; layer != 0; layer--) {
        Node &this_node = node_vec[this_node_id];
        Node &parent = node_vec[this_node.parent];
        size_t sister_id = (this_node_id & 1) ? (this_node_id + 1) : (this_node_id - 1);
        Node &sister = node_vec[sister_id];
        parent.leaf = this_node.val < sister.val ? this_node.leaf : sister.leaf;
        parent.val = this_node.val < sister.val ? this_node.val : sister.val;

        this_node_id = this_node.parent;
    }
    winner_id = root().leaf;
    return root().val;
}

template <typename T>
T
TournamentTree<T>::replace_run(size_t winner_id)
{
    return replace_run(max_val, winner_id);
}

template class TournamentTree<ibd_rec1_t>;
template class TournamentTree<uint16_t>;
