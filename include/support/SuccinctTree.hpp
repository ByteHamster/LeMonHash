#pragma once

#include "GrowingPastaBitVector.hpp"

/**
 * Balanced Parentheses Succinct Tree.
 * Allows to retrieve the rank of an inner nodes in order to compactly store data only for inner nodes.
 */
class SuccinctTree {
    private:
        static constexpr bool OPENING_NODE = false;
        static constexpr bool CLOSING_NODE = true;

        GrowingPastaBitVector succinctTreeRepresentation;
        GrowingPastaBitVector nodeIsInnerNode;
        pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES> *nodeIsInnerNodeRankSelect = nullptr;
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
            succinctTreeRepresentation.shrinkToFit();
            nodeIsInnerNodeRankSelect = new pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES>(nodeIsInnerNode.data);
        }

        size_t spaceBits() {
            return succinctTreeRepresentation.data.size()
                + nodeIsInnerNode.data.size()
                + 8 * nodeIsInnerNodeRankSelect->space_usage();
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
                    nextSibling();
                }
            }

            // Note that calling this on the last child renders the internal state invalid,
            // but the nodeId then still points to a DFS number one more than the size of the last child.
            void nextSibling() {
                size_t excess = 0;
                do {
                    if (tree.succinctTreeRepresentation.data[indexInTreeRepresentation] == OPENING_NODE) {
                        excess++;
                        nodeId++;
                    } else {
                        excess--;
                    }
                    indexInTreeRepresentation++;
                } while (excess != 0);
            }

            size_t innerNodeRank() {
                return tree.nodeIsInnerNodeRankSelect->rank1(nodeId);
            }
        };
};