#ifndef VARIANT_TYPES_H
#define VARIANT_TYPES_H

#include <string>
#include <functional> // std::hash のために必要

// フェージングの状態を表す列挙型
enum class Phase {
    Unphased,
    Hap1,
    Hap2
};

// ハッシュマップのキーとして使うための構造体
struct VariantKey {
    long long pos;
    std::string ref;
    std::string alt;

    bool operator==(const VariantKey& o) const {
        return pos == o.pos && ref == o.ref && alt == o.alt;
    }
};

// VariantKey を std::unordered_map で使うためのハッシュ関数
struct VariantKeyHash {
    std::size_t operator()(const VariantKey& k) const {
        std::size_t h1 = std::hash<long long>{}(k.pos);
        std::size_t h2 = std::hash<std::string>{}(k.ref);
        std::size_t h3 = std::hash<std::string>{}(k.alt);
        // ビットシフトしてXORをとる標準的なハッシュ合成
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

#endif // VARIANT_TYPES_H
