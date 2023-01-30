#pragma once

class TreePath {
    private:
        uint64_t hash = 88172645463325252ull;

        explicit TreePath(uint64_t hash) : hash(hash) {
        }

    public:
        TreePath() = default;

        [[nodiscard]] TreePath getChild(size_t childIndex) const {
            return TreePath(util::remix(hash + childIndex));
        }

        [[nodiscard]] uint64_t currentNodeHash() const {
            return hash;
        }

        [[nodiscard]] uint64_t alternativeHash() const {
            return util::remix(hash + 378368305212772073ul);
        }
};
