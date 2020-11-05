#ifndef ID_H
#define ID_H

#include <cstddef>
#include <functional>
#include <limits>

// Typesafe ID class such that there are compiler errors if different IDs are
// mixed. The template parameter T is just there to assure this behavior.
// Additionally, we have a member function which can check for validity.
template <typename T>
struct ID {
public:
    using IDType = std::size_t;
    static constexpr IDType invalid_value = std::numeric_limits<IDType>::max();

    ID(IDType nid = invalid_value) noexcept : id(nid) {}

    operator IDType() const {
        return id;
    }

    IDType operator+=(ID<T> other) {
        return id += other.id;
    }

    IDType operator-=(ID<T> other) {
        return id -= other.id;
    }

    IDType operator*=(ID<T> other) {
        return id *= other.id;
    }

    IDType operator++() {
        return ++id;
    }

    IDType operator--() {
        return --id;
    }

    bool operator!=(ID<T> other) const {
        return id != other.id;
    }

    bool valid() const {
        return id != invalid_value;
    }

    void invalidate() {
        id = invalid_value;
    }

private:
    IDType id;
};

// define custom hash function to be able to use IDs with maps/sets
namespace std
{
    template <typename T>
    struct hash<ID<T>>
    {
        using IDType = typename ID<T>::IDType;
        std::size_t operator()(ID<T> const& id) const noexcept
        {
            return std::hash<IDType>()(id);
        }
    };
} // std
#endif
