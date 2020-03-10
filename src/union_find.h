#pragma once

#include "id.h"
#include "random.h"

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>

template <typename T>
class UnionFind
{
public:
	using Value = T;
	using Values = typename std::vector<Value>;

	struct Element
	{
		ID<Element> const id;
		Value const value;

		ID<Element> parent_id;

		Element(ID<Element> id, Value value)
			: id(id), value(value) {}
	};
	using ElementID = ID<Element>;
	using Elements = std::vector<Element>;

	UnionFind(Values const& base_set);
	// returns the new root value
	Value uniteSets(Value value1, Value value2);
	// returns the root value of value
	Value findRoot(Value value);
	std::unordered_set<Value> const& getRoots() const;

private:
	Elements elements;
	std::unordered_set<Value> roots;

	std::unordered_map<Value, ElementID> to_id;
	Random random;
};

template <typename T>
UnionFind<T>::UnionFind(Values const& base_set)
{
	ElementID id = 0;
	for (auto value: base_set) {
		elements.emplace_back(id, value);
		roots.insert(value);
		to_id.emplace(value, id);

		++id;
	}
}

template <typename T>
auto UnionFind<T>::findRoot(Value value) -> Value
{
	auto id = to_id[value];

	// find root
	ElementID current_id = id;
	while (elements[current_id].parent_id.valid()) {
		current_id = elements[current_id].parent_id;
	}

	auto root = current_id;

	// path compression
	current_id = id;
	while (current_id != root) {
		auto parent_id = elements[current_id].parent_id;
		elements[current_id].parent_id = root;
		current_id = parent_id;
	}

	return elements[root].value;
}

template <typename T>
auto UnionFind<T>::uniteSets(Value value1, Value value2) -> Value
{
	auto root1_id = to_id[findRoot(value1)];
	auto root2_id = to_id[findRoot(value2)];
	auto& root1 = elements[root1_id];
	auto& root2 = elements[root2_id];

	// use randomization to decide which tree is hung under which, and delete
	// new child node from roots vector.
	if (random.throwCoin()) {
		root1.parent_id = root2_id;
		roots.erase(elements[root1_id].value);
		return elements[root2_id].value;
	}
	else {
		root2.parent_id = root1_id;
		roots.erase(elements[root2_id].value);
		return elements[root1_id].value;
	}
}

template <typename T>
auto UnionFind<T>::getRoots() const -> std::unordered_set<Value> const&
{
	return roots;
}
