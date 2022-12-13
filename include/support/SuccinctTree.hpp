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
        sdsl::bit_vector::rank_1_type succinctTreeRepresentationRank;
        sdsl::bp_support_g<> succinctTreeBpSupport;
        GrowingSdslBitVector nodeIsInnerNode;
        sdsl::bit_vector::rank_1_type nodeIsInnerNodeRank;
    public:
        SuccinctTree() = default;

        void openInnerNode() {
            succinctTreeRepresentation.append(OPENING_NODE);
            nodeIsInnerNode.append(true);
        }

        void closeInnerNode() {
            succinctTreeRepresentation.append(CLOSING_NODE);
        }

        void appendLeaf() {
            succinctTreeRepresentation.append(OPENING_NODE);
            succinctTreeRepresentation.append(CLOSING_NODE);
            nodeIsInnerNode.append(false);
        }

        void build() {
            nodeIsInnerNode.shrinkToFit();
            nodeIsInnerNodeRank = sdsl::bit_vector::rank_1_type(&nodeIsInnerNode.data);
            succinctTreeRepresentation.shrinkToFit();
            succinctTreeRepresentationRank = sdsl::bit_vector::rank_1_type(&succinctTreeRepresentation.data);
            succinctTreeBpSupport = sdsl::bp_support_g<>(&succinctTreeRepresentation.data);
        }

        size_t spaceBits() {
            return succinctTreeRepresentation.data.size()
                + nodeIsInnerNode.data.size()
                + nodeIsInnerNodeRank.bit_size()
                + succinctTreeBpSupport.bit_size()
                + succinctTreeRepresentationRank.bit_size();
        }

        struct Reader {
            size_t indexInTreeRepresentation = 0;
            size_t nodeId = 0;
            SuccinctTree &tree;

            explicit Reader(SuccinctTree &tree) : tree(tree) {
            }

            // TODO: This is possible in constant time. For now, just iterate.
            void skipToNthChild(size_t child) {
                indexInTreeRepresentation++;
                nodeId++;
                for (size_t i = 0; i < child; i++) {
                    indexInTreeRepresentation = tree.succinctTreeBpSupport.find_close(indexInTreeRepresentation) + 1;
                }
                nodeId = tree.succinctTreeRepresentationRank.rank(indexInTreeRepresentation);
            }

            // Note that calling this on the last child renders the internal state invalid,
            // but the nodeId then still points to a DFS number one more than the size of the last child.
            void nextSibling() {
                indexInTreeRepresentation = tree.succinctTreeBpSupport.find_close(indexInTreeRepresentation) + 1;
                nodeId = tree.succinctTreeRepresentationRank.rank(indexInTreeRepresentation);
            }

            size_t innerNodeRank() {
                return tree.nodeIsInnerNodeRank.rank(nodeId);
            }
        };
};