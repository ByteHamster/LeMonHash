#pragma once

#include "GrowingSdslBitVector.hpp"
#include "sdsl/bp_support.hpp"

/**
 * Balanced Parentheses Succinct Tree.
 * Allows to retrieve the rank of an inner nodes in order to compactly store data only for inner nodes.
 */
class SuccinctTree {
    private:
        static constexpr bool OPENING_NODE = true;
        static constexpr bool CLOSING_NODE = false;

        GrowingSdslBitVector succinctTreeRepresentation;
        sdsl::bp_support_gg<> succinctTreeBpSupport;
        sdsl::rank_support_v5<10, 2> leafNodeRank;
    public:
        SuccinctTree() = default;

        void openInnerNode() {
            succinctTreeRepresentation.append(OPENING_NODE);
        }

        void closeInnerNode() {
            succinctTreeRepresentation.append(CLOSING_NODE);
        }

        void appendLeaf() {
            succinctTreeRepresentation.append(OPENING_NODE);
            succinctTreeRepresentation.append(CLOSING_NODE);
        }

        void build() {
            succinctTreeRepresentation.shrinkToFit();
            succinctTreeBpSupport = sdsl::bp_support_gg<>(&succinctTreeRepresentation.data);
            leafNodeRank = sdsl::rank_support_v5<10, 2>(&succinctTreeRepresentation.data);
        }

        size_t spaceBits() {
            return succinctTreeRepresentation.data.size()
                   + leafNodeRank.bit_size()
                   + succinctTreeBpSupport.bit_size();
        }

        size_t spaceBitsWithoutIndices() {
            return succinctTreeRepresentation.data.size();
        }

        struct Reader {
            size_t indexInTreeRepresentation = 0;
            size_t nodeId = 0;
            SuccinctTree &tree;

            explicit Reader(SuccinctTree &tree) : tree(tree) {
            }

            void skipToNthChild(size_t child) {
                indexInTreeRepresentation++;
                nodeId++;
                for (size_t i = 0; i < child; i++) {
                    indexInTreeRepresentation = tree.succinctTreeBpSupport.find_close(indexInTreeRepresentation) + 1;
                }
                nodeId = tree.succinctTreeBpSupport.rank(indexInTreeRepresentation - 1);
            }

            // Note that calling this on the last child renders the internal state invalid,
            // but the nodeId then still points to a DFS number one more than the size of the last child.
            void nextSibling() {
                indexInTreeRepresentation = tree.succinctTreeBpSupport.find_close(indexInTreeRepresentation) + 1;
                nodeId = tree.succinctTreeBpSupport.rank(indexInTreeRepresentation - 1);
            }

            size_t innerNodeRank() {
                return nodeId - tree.leafNodeRank.rank(indexInTreeRepresentation);
            }
        };
};