#include <algorithm>
#include <cassert>
#include <complex>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

const int BASE = 1e9;
const int WIDTH = 9;
const long double PI = std::acos(-1.0L);

typedef std::complex<long double> complex;

class BigInt {
public:
    std::vector<int> digits_;

    BigInt (int64_t number = 0) {
        assert(number >= 0);
        do {
            digits_.push_back(static_cast<int>(number % BASE));
            number /= BASE;
        } while (number > 0);
        RemoveLeadingZeros();
    }

    BigInt(const std::string &str) {
        const int size = static_cast<int>(str.size());
        for (int index_of_group = 1, count_of_groups = size / WIDTH; index_of_group <= count_of_groups; ++index_of_group) {
            digits_.push_back(std::stoi(str.substr(size - index_of_group * WIDTH, WIDTH)));
        }
        if (size % WIDTH != 0) {
            digits_.push_back(std::stoi(str.substr(0, size % WIDTH)));
        }
        RemoveLeadingZeros();
    }

    BigInt(const std::vector<int> &digits) : digits_(digits) {
        RemoveLeadingZeros();
    }

    BigInt& RemoveLeadingZeros() {
        while (digits_.back() == 0 && static_cast<int>(digits_.size()) > 1) digits_.pop_back();
        for (auto digit: digits_) assert(0 <= digit && digit < BASE);
        return *this;
    }

    int Compare(const BigInt& other) const;
    BigInt SlowMultiplication(const BigInt& other) const;
    BigInt FastMultiplication(const BigInt& other) const;
    BigInt Multiplication(const BigInt& other) const;
    std::pair<BigInt, BigInt> DivideMod(const BigInt& other) const;

    friend std::istream& operator>>(std::istream& in, BigInt& number) {
        std::string str;
        in >> str;
        number = BigInt(str);
        return in;
    }

    friend std::ostream& operator<<(std::ostream& out, const BigInt& number) {
        out << number.digits_.back();
        for (int i = static_cast<int>(number.digits_.size()) - 2; i >= 0; --i) {
            out << std::setw(WIDTH) << std::setfill('0') << number.digits_[i];
        }
        return out << std::setfill(' ');
    }

    std::string ToString() {
        std::stringstream ss;
        ss << *this;
        return ss.str();
    }

    BigInt& operator+=(const int number) {
        assert(number >= 0);
        if (number >= BASE) {
            return *this += BigInt(number);
        }
        int remainder = number;
        for (int i = 0; remainder > 0; ++i) {
            if (i >= static_cast<int>(digits_.size())) digits_.push_back(0);
            remainder += digits_[i];
            if (remainder >= BASE) {
                digits_[i] = remainder - BASE;
                remainder = 1;
            } else {
                digits_[i] = remainder;
                remainder = 0;
            }
        }
        return this->RemoveLeadingZeros();
    }

    BigInt& operator+=(const BigInt& other) {
        if (other.digits_.size() == 1u) {
            return *this += other.digits_[0];
        }
        const int size1 = static_cast<int>(this->digits_.size());
        const int size2 = static_cast<int>(other.digits_.size());
        int remainder = 0;
        for (int i = 0; i < std::max(size1, size2) || remainder > 0; ++i) {
            int div1 = i < size1 ? this->digits_[i] : (digits_.push_back(0), 0);
            int div2 = i < size2 ? other.digits_[i] : 0;
            remainder += div1 + div2;
            auto divisor = remainder / BASE;
            digits_[i] = remainder - divisor * BASE;
            remainder = divisor;
        }
        return this->RemoveLeadingZeros();
    }


    BigInt& operator-=(const int number) {
        assert(number >= 0);
        if (number >= BASE) {
            return *this -= BigInt(number);
        }
        int remainder = -number;
        for (int i = 0; i < (int) digits_.size() && remainder < 0; ++i) {
            remainder += digits_[i];
            if (remainder < 0) {
                digits_[i] = remainder + BASE;
                remainder = -1;
            } else {
                digits_[i] = remainder;
                remainder = 0;
            }
        }
        assert(remainder == 0);
        return this->RemoveLeadingZeros();
    }

    BigInt& operator-=(const BigInt& other) {
        if (other.digits_.size() == 1u) {
            return *this -= other.digits_[0];
        }
        const int s1 = static_cast<int>(this->digits_.size());
        const int s2 = static_cast<int>(other.digits_.size());
        assert(s1 >= s2);
        int remainder = 0;
        for (int i = 0; i < s1; ++i) {
            int d2 = i < s2 ? other.digits_[i] : 0;
            remainder += this->digits_[i] - d2;
            if (remainder < 0) {
                digits_[i] = remainder + BASE;
                remainder = -1;
            } else {
                digits_[i] = remainder;
                remainder = 0;
                if (i >= s2) break;
            }
        }
        assert(remainder == 0);
        return this->RemoveLeadingZeros();
    }

    BigInt& operator*=(const unsigned int number) {
        assert(number >= 0);
        if (number >= BASE) {
            return *this *= BigInt(number);
        }
        int64_t remainder = 0;
        for (auto &digit: digits_) {
            remainder += 1LL * digit * number;
            auto divisor = remainder / BASE;
            digit = remainder - divisor * BASE;
            remainder = divisor;
        }
        if (remainder > 0) digits_.push_back(remainder);
        return this->RemoveLeadingZeros();
    }

    BigInt& operator*=(const BigInt& other);
    BigInt& operator/=(const int num);
    BigInt& operator/=(const BigInt& other);
    BigInt& operator%=(const BigInt& other);
};

BigInt operator+(const BigInt&, const BigInt&);
BigInt operator-(const BigInt&, const BigInt&);
BigInt operator*(const BigInt&, const BigInt&);
BigInt operator/(const BigInt&, const BigInt&);
BigInt operator%(const BigInt&, const BigInt&);

BigInt operator+(const BigInt&, const int);
BigInt operator+(const int, const BigInt&);
BigInt operator-(const BigInt&, const int);
BigInt operator*(const BigInt&, const int);
BigInt operator*(const int, const BigInt&);
BigInt operator/(const BigInt&, const int);

bool operator<(const BigInt&, const BigInt&);
bool operator>(const BigInt&, const BigInt&);
bool operator<=(const BigInt&, const BigInt&);
bool operator>=(const BigInt&, const BigInt&);
bool operator==(const BigInt&, const BigInt&);
bool operator!=(const BigInt&, const BigInt&);

BigInt BigInt::SlowMultiplication(const BigInt& other) const {
    if (other.digits_.size() == 1u) {
        return *this * other.digits_[0];
    }
    const int size1 = static_cast<int>(this->digits_.size());
    const int size2 = static_cast<int>(other.digits_.size());
    std::vector<int> temporary(size1 + size2);
    for (int i = 0; i < size1; ++i) {
        int64_t remainder = 0;
        for (int j = 0; j < size2; ++j) {
            remainder += temporary[i + j] + 1LL * this->digits_[i] * other.digits_[j];
            auto divisor = remainder / BASE;
            temporary[i + j] = remainder - divisor * BASE;
            remainder = divisor;
        }
        if (remainder > 0) temporary[i + size2] += remainder;
    }
    return BigInt(temporary);
}

BigInt BigInt::FastMultiplication(const BigInt& other) const {
    if (other.digits_.size() == 1u) {
        return *this * other.digits_[0];
    }

    std::function<int(int, int)> reverse = [](int number, int count_of_bits) {
        int result = 0;
        for (int i = 0; i < count_of_bits; ++i) {
            if (number & (1 << i)) {
                result |= 1 << (count_of_bits - 1 - i);
            }
        }
        return result;
    };

    std::function<void(std::vector<complex>&, bool)> fft = [&reverse](std::vector<complex>& a, bool invert) {
        const int n = static_cast<int>(a.size());
        int count_of_bits = 0;
        while ((1 << count_of_bits) < n) ++count_of_bits;

        for (int i = 0; i < n; ++i) {
            if (i < reverse(i, count_of_bits)) {
                std::swap(a[i], a[reverse(i, count_of_bits)]);
            }
        }

        for (int length = 2; length <= n; length <<= 1) {
            auto ang = 2 * PI / length * (invert ? -1 : 1);
            complex wlen(std::cos(ang), std::sin(ang));
            for (int i = 0; i < n; i += length) {
                complex w(1);
                for (int j = 0; j < length / 2; ++j) {
                    complex u = a[i + j];
                    complex v = a[i + j + length / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + length / 2] = u - v;
                    w *= wlen;
                }
            }
        }
        if (invert) {
            for (int i = 0; i < n; ++i) {
                a[i] /= n;
            }
        }
    };

    assert(BASE == 1000 * 1000 * 1000);
    std::function<std::vector<complex>(const BigInt&)> prepare = [](const BigInt& number) {
        std::vector<complex> result;
        result.reserve(3 * number.digits_.size());
        for (auto digit: number.digits_) {
            result.push_back(digit % 1000);
            result.push_back(digit / 1000 % 1000);
            result.push_back(digit / 1000000);
        }
        return result;
    };

    auto fa = prepare(*this);
    auto fb = prepare(other);

    int n = 1;
    while (n < static_cast<int>(std::max(fa.size(), fb.size()))) {
        n *= 2;
    }
    n *= 2;
    fa.resize(n);
    fb.resize(n);

    fft(fa, false);
    fft(fb, false);

    for (int i = 0; i < n; ++i) {
        fa[i] *= fb[i];
    }

    fft(fa, true);

    std::vector<int64_t> temporary(n);
    for (int i = 0; i < static_cast<int>(fa.size()); ++i) {
        temporary[i] = static_cast<int64_t>(fa[i].real() + 0.5);
    }

    int64_t carry = 0;
    for (int i = 0; i < n || carry > 0; ++i) {
        if (i >= n) temporary.push_back(0);
        temporary[i] += carry;
        carry = temporary[i] / 1000;
        temporary[i] -= carry * 1000;
        assert(temporary[i] >= 0);
    }

    std::vector<int> result;
    result.reserve(this->digits_.size() + other.digits_.size());

    for (int i = 0; i < n; i += 3) {
        int c = temporary[i];
        int b = i + 1 < n ? temporary[i + 1] : 0;
        int a = i + 2 < n ? temporary[i + 2] : 0;
        result.push_back(c + 1000 * (b + 1000 * a));
    }
    return BigInt(result);
}

BigInt BigInt::Multiplication(const BigInt& other) const {
    int length1 = static_cast<int>(this->digits_.size());
    int length2 = static_cast<int>(other.digits_.size());
    int temporary = 3 * std::max(length1, length2);
    int power = 1;
    while (power < temporary) power *= 2;
    power *= 2;
    int op1 = length1 * length2;
    int op2 = 3 * power * std::log(power) / std::log(2);
    return op1 >= 15 * op2 ? FastMultiplication(other) : SlowMultiplication(other);
}

BigInt& BigInt::operator/=(const int number) {
    assert(number > 0);
    if (number >= BASE) {
        return *this /= BigInt(number);
    }
    int64_t remainder = 0;
    for (int j = (int) digits_.size() - 1; j >= 0; --j) {
        (remainder *= BASE) += digits_[j];
        auto divisor = remainder / number;
        digits_[j] = divisor;
        remainder -= divisor * number;
    }
    return this->RemoveLeadingZeros();
}


int operator%(const BigInt& a, const unsigned int number) {
    assert(number > 0);
    int64_t remainder = 0;
    for (int i = static_cast<int>(a.digits_.size()) - 1; i >= 0; --i) {
        ((remainder *= BASE) += a.digits_[i]) %= number;
    }
    return remainder;
}

std::pair<BigInt, BigInt> BigInt::DivideMod(const BigInt& other) const {
    if (other.digits_.size() == 1u) {
        return {std::move(*this / other.digits_[0]), *this % other.digits_[0]};
    }
    const int norm = BASE / (other.digits_.back() + 1);
    const BigInt a = *this * norm;
    const BigInt b = other * norm;
    const int a_size = static_cast<int>(a.digits_.size());
    const int b_size = static_cast<int>(b.digits_.size());
    BigInt q, r;
    q.digits_.resize(a_size);
    for (int i = a_size - 1; i >= 0; --i) {
        r *= BASE;
        r += a.digits_[i];
        int s1 = static_cast<int>(r.digits_.size()) <= b_size ? 0 : r.digits_[b_size];
        int s2 = static_cast<int>(r.digits_.size()) <= b_size - 1 ? 0 : r.digits_[b_size - 1];
        int d = (1LL * BASE * s1 + s2) / b.digits_.back();
        auto temp = b * d;
        while (r < temp) {
            r += b;
            --d;
        }
        r -= temp;
        q.digits_[i] = d;
    }
    return {std::move(q.RemoveLeadingZeros()), std::move(r /= norm)};
}

int BigInt::Compare(const BigInt& other) const {
    if (this->digits_.size() > other.digits_.size()){
        return 1;
    }
    if (this->digits_.size() < other.digits_.size()) {
        return -1;
    }
    for (int i = (int) digits_.size() - 1; i >= 0; --i) {
        if (this->digits_[i] > other.digits_[i]) {
            return 1;
        }
        if (this->digits_[i] < other.digits_[i]) {
            return -1;
        }
    }
    return 0;
}

bool operator<(const BigInt& a, const BigInt& b) {
    return a.Compare(b) < 0;
}

bool operator>(const BigInt& a, const BigInt& b) {
    return a.Compare(b) > 0;
}

bool operator==(const BigInt& a, const BigInt& b) {
    return a.Compare(b) == 0;
}

bool operator<=(const BigInt& a, const BigInt& b) {
    return a.Compare(b) <= 0;
}

bool operator>=(const BigInt& a, const BigInt& b) {
    return a.Compare(b) >= 0;
}

bool operator!=(const BigInt& a, const BigInt& b) {
    return a.Compare(b) != 0;
}

BigInt operator+(const BigInt& a, const BigInt& b) {
    return BigInt(a) += b;
}

BigInt operator-(const BigInt& a, const BigInt& b) {
    return BigInt(a) -= b;
}

BigInt operator*(const BigInt& a, const BigInt& b) {
    return a.Multiplication(b);
}

BigInt operator/(const BigInt& a, const BigInt& b) {
    return a.DivideMod(b).first;
}

BigInt operator%(const BigInt& a, const BigInt& b) {
    return a.DivideMod(b).second;
}

BigInt& BigInt::operator*=(const BigInt& other) {
    return *this = *this * other;
}

BigInt& BigInt::operator/=(const BigInt& other) {
    return *this = *this / other;
}

BigInt& BigInt::operator%=(const BigInt& other) {
    return *this = *this % other;
}

BigInt operator+(const BigInt& a, const int b) {
    return BigInt(a) += b;
}

BigInt operator+(const int a, const BigInt& b) {
    return b + a;
}

BigInt operator-(const BigInt& a, const int b) {
    return BigInt(a) -= b;
}

BigInt operator*(const BigInt& a, const int b) {
    return BigInt(a) *= b;
}

BigInt operator*(const int a, const BigInt& b) {
    return b * a;
}

BigInt operator/(const BigInt& a, const int b) {
    return BigInt(a) /= b;
}

BigInt Power(BigInt a, BigInt n, BigInt p) {
    BigInt res = 1;
    while (n > 0) {
        if (n % 2 != 0) res *= a;
        a *= a;
        a %= p;
        n /= 2;
        res %= p;
    }
    return res;
}

char IntToChar(int value) {
    if (value >= 0 && value <= 9) {
        return static_cast<char>(value + 48);
    }
    if (value >= 10 && value <= 35) {
        return static_cast<char>(value + 55);
    }
    if (value >= 36 && value <= 61) {
        return static_cast<char>(value + 61);
    }
    if (value == 62) {
        return ' ';
    }
    if (value == 63) {
        return '.';
    }
    return '$';
}

int main() {
    BigInt p, k;
    std::cin >> p >> k;
    BigInt number = 0;
    BigInt power = 1;
    int64_t X, Y;
    while (std::cin >> X >> Y) {
        BigInt x(X);
        BigInt y(Y);
        BigInt result = y * (Power(x, p - k - 1, p)) % p;
        number += result * power;
        power *= p;
    }
    while (number != 0) {
        std::cout << IntToChar(number % 64);
        number /= 64;
    }
    
    return 0;
}
