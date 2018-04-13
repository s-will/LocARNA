#ifndef LOCARNA_NAMED_ARGUMENTS_HH
#define LOCARNA_NAMED_ARGUMENTS_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>

/**
 * Mechanism to name arguments via unique argument classes;
 * in the function, the actual arguments can be retrieved from an argument pack
 * by their type.
 *
 * Requires functions with named arguments to use (variadic) templates
 */

#include <tuple>
#include <type_traits>

namespace LocARNA {

    // template meta program to test whether a tuple has a type
    template <typename T, typename Tuple>
    struct has_type;

    template <typename T>
    struct has_type<T, std::tuple<>> : std::false_type {};

    template <typename T, typename U, typename... Ts>
    struct has_type<T, std::tuple<U, Ts...>> : has_type<T, std::tuple<Ts...>> {
    };

    template <typename T, typename... Ts>
    struct has_type<T, std::tuple<T, Ts...>> : std::true_type {};

    // check whether T is contained in Us...
    template <typename T, typename... Us>
    struct contains;

    template <typename T>
    struct contains<T, std::tuple<>> : std::false_type {};

    template <typename T, typename U, typename... Us>
    struct contains<T, std::tuple<U, Us...>> : contains<T, std::tuple<Us...>> {
    };

    template <typename T, typename... Us>
    struct contains<T, std::tuple<T, Us...>> : std::true_type {};

    template <typename T, typename... Us>
    using contains_t = typename contains<T, Us...>::type;

    // check that types in Ts are a subset of Us...
    template <typename Tuple1, typename Tuple2>
    struct type_subset_of;

    template <typename... Ts>
    struct type_subset_of<std::tuple<Ts...>, std::tuple<>> : std::false_type {};

    template <typename... Us>
    struct type_subset_of<std::tuple<>, std::tuple<Us...>> : std::true_type {};

    template <typename T, typename... Ts, typename... Us>
    struct type_subset_of<std::tuple<T, Ts...>, std::tuple<Us...>> {
        static constexpr bool value =
            type_subset_of<std::tuple<Ts...>, std::tuple<Us...>>::value &&
            contains<T, std::tuple<Us...>>::value;
    };

    // get with default value
    //
    template <
        typename T,
        typename... Args,
        std::enable_if_t<!has_type<T, std::tuple<Args...>>::value> * = nullptr>
    constexpr T
    get_def(const std::tuple<Args...> &t, const T &def) noexcept {
        return def;
    }

    template <
        typename T,
        typename... Args,
        std::enable_if_t<has_type<T, std::tuple<Args...>>::value> * = nullptr>
    constexpr const T &
    get_def(const std::tuple<Args...> &t, const T &def) noexcept {
        return std::get<T>(t);
    }

    template <typename T>
    class NamedArgument {
    public:
        using value_type = T;
        NamedArgument() : value_(){};
        NamedArgument(const value_type &value) : value_(value){};
        NamedArgument(value_type &&value) : value_(value){};

        value_type
        value() const {
            return value_;
        }

    private:
        value_type value_;
    };

    //! define a named argument with name and type
#   define DEFINE_NAMED_ARG(name, type)                     \
    struct name : public NamedArgument<type> {           \
        name(type value) : NamedArgument<type>(value) {} \
    }

    //! define a named argument with name, type, and default value
#   define DEFINE_NAMED_ARG_DEFAULT(name, type, default_value) \
    struct name : public NamedArgument<type> {              \
        name() : NamedArgument<type>(default_value) {}      \
        name(type value) : NamedArgument<type>(value) {}    \
    }

    //! @brief get a mandatory named argument from argument pack;
    //! ensure that argument is given
    template <typename name, typename ArgTuple>
    auto
    get_named_arg(const ArgTuple &args) {
        static_assert(has_type<name, ArgTuple>::value,
                      "Required named argument missing.");
        return std::get<name>(args).value();
    }

    // //! get a named argument from argument pack; return default if argument does not
    // //! exist
    // template <typename name, typename ValueType, typename ArgsTuple>
    // auto
    // get_named_arg_def(ValueType &&def, const ArgsTuple &args) {
    //     return get_def<name>(args,
    //                          name(std::forward<ValueType>(def))).value();
    // }

    //! get an optional named argument from argument pack; return default value
    //! (empty
    //! constructor of named parameter class) if argument does not
    //! exist
    template <typename name, typename ArgsTuple>
    auto
    get_named_arg_opt(const ArgsTuple &args) {
        return get_def<name>(args,
                             name()).value();
    }

    //! define a named argument with name and type; additionally declare variable
    //! (to be used as class feature)
#   define DEFINE_NAMED_ARG_FEATURE(name, type)  \
    type name##_;                                \
    DEFINE_NAMED_ARG(name, type)

    //! define a named argument with name, type, and default; additionally declare variable
    //! (to be used as class feature)
#   define DEFINE_NAMED_ARG_DEFAULT_FEATURE(name, type, default_value) \
    type name##_;                                                   \
    DEFINE_NAMED_ARG_DEFAULT(name, type, default_value)

} // end namespace LocARNA

#endif // LOCARNA_NAMED_ARGUMENTS_HH
