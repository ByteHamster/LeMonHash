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

        void toDotFile() {
            size_t maxExpectedDepth = 35;
            float scaleY = 200;
            float scaleX = (maxExpectedDepth * scaleY) / (succinctTreeRepresentation.size / 2); // Try making it square
            bool centerParentOnChildren = true;

            std::cout<<"digraph {"<<std::endl;
            std::vector<size_t> stack;
            std::vector<size_t> childStack;
            stack.push_back(-1);
            childStack.push_back(0);
            size_t dfsNum = 0;
            for (size_t i = 0; i < succinctTreeRepresentation.data.size(); i++) {
                if (succinctTreeRepresentation.data[i] == OPENING_NODE) {
                    dfsNum++;
                    if (i != 0) {
                        std::cout<<"  "<<stack.back()<<" -> "<<dfsNum<<std::endl;
                    }
                    stack.push_back(dfsNum);
                    childStack.back()++;
                    childStack.push_back(0);
                } else {
                    size_t children = childStack.back();
                    size_t endDfsNum = stack.back();
                    std::cout<<"  "<<endDfsNum<<" [ "<<std::endl;
                    std::cout<<"    pos = \""<<scaleX * (endDfsNum + (centerParentOnChildren ? children/2 : 0))
                                        <<","<<scaleY * stack.size()<<"\""<<std::endl;
                    std::cout<<"  ]"<<std::endl;
                    stack.pop_back();
                    childStack.pop_back();
                    childStack.back() += children;
                }
            }
            std::cout<<"}"<<std::endl;
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