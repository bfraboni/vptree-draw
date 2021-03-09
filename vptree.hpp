// Based on http://stevehanov.ca/blog/index.php?id=130 by Steve Hanov

#pragma once

#include <algorithm>
#include <vector>
#include <queue>
#include <limits>
#include <random>
#include <cmath>
#include <stdexcept>

namespace vpt {

typedef std::vector<double> Vector;
typedef std::function<double(const Vector& v1, const Vector& v2)> Metric;
typedef std::pair<std::vector<double>, std::vector<int>> DistancesIndices;
typedef std::pair<std::vector<std::vector<double>>, std::vector<std::vector<int>>> BatchDistancesIndices;


template<class InputIterator>
double sum(InputIterator begin, InputIterator end) {
    double result = 0;
    for (; begin != end; ++begin) {
        result += *begin;
    }
    return result;
}

template<typename Container>
struct EuclideanMetric {
    double operator() (const Container& v1, const Container& v2) const {
        Vector diffSquares(v1.size());
        std::transform(v1.begin(), v1.end(), v2.begin(), diffSquares.begin(), [] (double lhs, double rhs) {
            return (lhs - rhs) * (lhs - rhs);
        });
        auto sum = vpt::sum(diffSquares.begin(), diffSquares.end());
        return std::sqrt(sum);
    }
};

class DimensionMismatch: public std::runtime_error {
public:
    DimensionMismatch(int expected, int got)
    : std::runtime_error("Item dimension doesn't match: expected " + std::to_string(expected) + ", got " + std::to_string(got))
    {}
};

class Searcher;


class VpTree {
public:
    template<typename InputIterator>
    explicit VpTree(InputIterator start, InputIterator end, Metric metric = EuclideanMetric<Vector>());

    template<typename Container>
    explicit VpTree(const Container& container, Metric metric = EuclideanMetric<Vector>());
    explicit VpTree(std::initializer_list<Vector> list, Metric metric = EuclideanMetric<Vector>());

    DistancesIndices getNearestNeighbors(const Vector& target, int neighborsCount) const;
    template<typename VectorLike>
    DistancesIndices getNearestNeighbors(const VectorLike& target, int neighborsCount) const;
    DistancesIndices getNearestNeighbors(std::initializer_list<double> target, int neighborsCount) const;

    template<typename Container>
    BatchDistancesIndices getNearestNeighborsBatch(const Container& targets, int neighborsCount) const;
    BatchDistancesIndices getNearestNeighborsBatch(std::initializer_list<Vector> targets, int neighborsCount) const;

    const Metric getDistance;

public:
    struct Node {
        static const int Leaf = -1;

        Node(int item, double threshold = 0., int left = Leaf, int right = Leaf)
        : item(item), threshold(threshold), left(left), right(right)
        { }

        int item;
        double threshold;
        int left;
        int right;
    };

public:
    typedef std::pair<Vector, int> Item;

    std::vector<Item> items_;
    std::vector<Node> nodes_;

    std::mt19937 rng_;

    int dimension_;

public:
    template<typename InputIterator>
    std::vector<Item> makeItems(InputIterator start, InputIterator end);

    int makeTree(int lower, int upper);
    void selectRoot(int lower, int upper);
    void partitionByDistance(int lower, int pos, int upper);
    int makeNode(int item);
    Node root() const { return nodes_[0]; }

    friend class Searcher;
};

class Searcher {
public:
    typedef typename VpTree::Node Node;

public:
    explicit Searcher(const VpTree* tree, const Vector& target, int neighborsCount);

    std::pair<std::vector<double>, std::vector<int>> search();

    struct HeapItem {
        bool operator < (const HeapItem& other) const {
            return dist < other.dist;
        }

        int item;
        double dist;
    };

public:
    void searchInNode(const Node& node);

    const VpTree* tree_;
    Vector target_;
    int neighborsCount_;
    double tau_;
    std::priority_queue<HeapItem> heap_;
};



template<typename InputIterator>
VpTree::VpTree(InputIterator start, InputIterator end, Metric metric)
: getDistance(metric), items_(makeItems(start, end)), nodes_(), rng_() {
    std::random_device rd;
    rng_.seed(rd());
    nodes_.reserve(items_.size());
    makeTree(0, items_.size());
}

template<typename Container>
VpTree::VpTree(const Container& container, Metric metric)
: VpTree(container.begin(), container.end(), metric)
{ }

VpTree::VpTree(std::initializer_list<Vector> list, Metric metric)
: VpTree(list.begin(), list.end(), metric)
{ }

int VpTree::makeTree(int lower, int upper) {
    if (lower >= upper) {
        return Node::Leaf;
    } else if (lower + 1 == upper) {
        return makeNode(lower);
    } else {
        selectRoot(lower, upper);
        int median = (upper + lower) / 2;
        partitionByDistance(lower, median, upper);
        auto node = makeNode(lower);
        nodes_[node].threshold = getDistance(items_[lower].first, items_[median].first);
        nodes_[node].left = makeTree(lower + 1, median);
        nodes_[node].right = makeTree(median, upper);
        return node;
    }
}

void VpTree::selectRoot(int lower, int upper) {
    std::uniform_int_distribution<int> uni(lower, upper - 1);
    int root = uni(rng_);
    std::swap(items_[lower], items_[root]);
}

void VpTree::partitionByDistance(int lower, int pos, int upper) {
    std::nth_element(
        items_.begin() + lower + 1,
        items_.begin() + pos,
        items_.begin() + upper,
        [lower, this] (const Item& i1, const Item& i2) {
            return getDistance(items_[lower].first, i1.first) < getDistance(items_[lower].first, i2.first);
        });
}

int VpTree::makeNode(int item) {
    nodes_.push_back(Node(item));
    return nodes_.size() - 1;
}

template<typename InputIterator>
std::vector<std::pair<Vector, int>> VpTree::makeItems(InputIterator begin, InputIterator end) {
    if (begin != end) {
        dimension_ = begin->size();
    } else {
        dimension_ = -1;
    }

    std::vector<Item> res;
    for (int i = 0; begin != end; ++begin, ++i) {
        auto vec = Vector(begin->begin(), begin->end());
        res.push_back(std::make_pair(vec, i));

        auto lastDimension = res.back().first.size();
        if (lastDimension != dimension_) {
            throw DimensionMismatch(dimension_, lastDimension);
        }
    }
    return res;
}

template<typename VectorLike>
DistancesIndices VpTree::getNearestNeighbors(const VectorLike& target, int neighborsCount) const {
    return getNearestNeighbors(Vector(target.begin(), target.end()), neighborsCount);
}

DistancesIndices VpTree::getNearestNeighbors(std::initializer_list<double> target, int neighborsCount) const {
    return getNearestNeighbors(Vector(target.begin(), target.end()), neighborsCount);
}

DistancesIndices VpTree::getNearestNeighbors(const Vector& target, int neighborsCount) const {
    auto targetDimension = target.size();
    if (targetDimension != dimension_) {
        throw DimensionMismatch(dimension_, targetDimension);
    }
    Searcher searcher(this, target, neighborsCount);
    return searcher.search();
}

template<typename Container>
BatchDistancesIndices VpTree::getNearestNeighborsBatch(const Container& targets, int neighborsCount) const {
    std::vector<std::vector<double>> batchDistances(targets.size());
    std::vector<std::vector<int>> batchIndices(targets.size());
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < targets.size(); ++i) {
        std::tie(batchDistances[i], batchIndices[i]) = getNearestNeighbors(targets[i], neighborsCount);
    }
    return BatchDistancesIndices(batchDistances, batchIndices);
}

BatchDistancesIndices VpTree::getNearestNeighborsBatch(std::initializer_list<Vector> targets, int neighborsCount) const {
    return getNearestNeighborsBatch(std::vector<Vector>(targets.begin(), targets.end()), neighborsCount);
}


Searcher::Searcher(const VpTree* tree, const Vector& target, int neighborsCount)
: tree_(tree), target_(target), neighborsCount_(neighborsCount), tau_(std::numeric_limits<double>::max()), heap_()
{ }

DistancesIndices Searcher::search() {
    searchInNode(tree_->root());

    DistancesIndices results;
    while(!heap_.empty()) {
        results.first.push_back(heap_.top().dist);
        results.second.push_back(tree_->items_[heap_.top().item].second);
        heap_.pop();
    }
    std::reverse(results.first.begin(), results.first.end());
    std::reverse(results.second.begin(), results.second.end());
    return results;
}

void Searcher::searchInNode(const Node& node) {
    double dist = tree_->getDistance(tree_->items_[node.item].first, target_);

    if (dist < tau_) {
        if (heap_.size() == neighborsCount_)
            heap_.pop();

        heap_.push(HeapItem{node.item, dist});

        if (heap_.size() == neighborsCount_)
            tau_ = heap_.top().dist;
    }

    if (dist < node.threshold) {
        if (node.left != Node::Leaf && dist - tau_ <= node.threshold)
            searchInNode(tree_->nodes_[node.left]);

        if (node.right != Node::Leaf && dist + tau_ >= node.threshold)
            searchInNode(tree_->nodes_[node.right]);
    } else {
        if (node.right != Node::Leaf && dist + tau_ >= node.threshold)
            searchInNode(tree_->nodes_[node.right]);

        if (node.left != Node::Leaf && dist - tau_ <= node.threshold)
            searchInNode(tree_->nodes_[node.left]);
    }
}

}