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

BigInt p;

BigInt Power(BigInt a, BigInt n, BigInt p) {
    BigInt result = 1;
    while (n > 0) {
        if (n % 2 != 0) result *= a;
        a *= a;
        a %= p;
        n /= 2;
        result %= p;
    }
    return result;
}

BigInt Power(BigInt a, BigInt n) {
    BigInt result = 1;
    while (n > 0) {
        if (n % 2 != 0) result *= a;
        a *= a;
        a %= p;
        n /= 2;
        result %= p;
    }
    return result;
}

class Poly {
public:
    Poly(std::vector<BigInt> coefficients) : coefficients_(coefficients) {
        DeleteLeadingZeros();
    }

    Poly(const BigInt& coefficient = BigInt()) {
        coefficients_.push_back(coefficient);
        DeleteLeadingZeros();
    }

    size_t Degree() const {
        return coefficients_.size();
    }

    bool operator==(const Poly& other) const {
        return coefficients_ == other.coefficients_;
    }

    bool operator!=(const Poly& other) const {
        return coefficients_ != other.coefficients_;
    }

    Poly &operator+=(const Poly& other) {
        coefficients_.resize(std::max(coefficients_.size(), other.coefficients_.size()));
        for (int64_t i = 0; i != static_cast<int64_t>(other.coefficients_.size()); ++i) {
            coefficients_[i] += other.coefficients_[i];
        }
        DeleteLeadingZeros();
        return *this;
    }

    Poly operator+(const Poly& other) const {
        Poly result = *this;
        result += other;
        return result;
    }

    Poly operator-(Poly& other) const {
        Poly result = *this;
        result.DeleteLeadingZeros();
        other.DeleteLeadingZeros();
        for (int64_t i = 0; i != static_cast<int64_t>(other.coefficients_.size()); ++i) {
            result.coefficients_[i] = (result.coefficients_[i] + p - other.coefficients_[i]) % p;
        }
        result.DeleteLeadingZeros();
        return result;
    }

    const BigInt operator[](size_t i) const {
        if (i >= coefficients_.size()) {
            return 0;
        }
        return coefficients_[i];
    }

    template<typename Iterator>
    Poly(Iterator begin, Iterator end) : coefficients_(begin, end) {
        DeleteLeadingZeros();
    }

    typename std::vector<BigInt>::const_iterator begin() const {
        return coefficients_.cbegin();
    }

    typename std::vector<BigInt>::const_iterator end() const {
        return coefficients_.cend();
    }

    Poly &operator*=(const Poly& other) {
        std::vector<BigInt> result(coefficients_.size() + other.coefficients_.size() - 1);
        for (int64_t i = 0; i != static_cast<int64_t>(coefficients_.size()); ++i) {
            for (int64_t j = 0; j != static_cast<int64_t>(other.coefficients_.size()); ++j) {
                result[i + j] = (result[i + j]  + (coefficients_[i] * other.coefficients_[j]) % p) % p;
            }
        }
        coefficients_ = result;
        return *this;
    }

    Poly operator*(const Poly& other) {
        Poly result = *this;
        result *= other;
        return result;
    }

    Poly operator/(const Poly& other) {
        Poly polynomial1 = *this;
        Poly polynomial2;
        Poly result;
        std::vector<BigInt> coefficients(polynomial1.coefficients_.size(), 0);
        if (polynomial1.Degree() < other.Degree()) {
            result = Poly(BigInt(1));
            return result;
        }
        BigInt multiplier1;
        BigInt multiplier2;
        Poly polynomial3;
        while (polynomial1.Degree() >= other.Degree()) {
            polynomial2 = other;
            multiplier1 = Power(polynomial2.coefficients_.back(), p - 2);
            multiplier2 = (polynomial1.coefficients_.back() * multiplier1) % p;
            coefficients[polynomial1.Degree() - polynomial2.Degree()] = multiplier2;
            polynomial3 = (polynomial2 * Poly(multiplier2));
            std::vector<BigInt> div(polynomial1.Degree() - other.Degree() + 1);
            div.back() = 1;
            polynomial3 = polynomial3 * Poly(div);
            polynomial1 = polynomial1 - polynomial3;
        }
        result.coefficients_ = coefficients;
        result.DeleteLeadingZeros();
        return result;
    }

    Poly operator%(const Poly& other) {
        Poly polynomial1 = *this;
        polynomial1.DeleteLeadingZeros();
        Poly result;
        if (polynomial1.Degree() < other.Degree()) {
            return polynomial1;
        }
        Poly polynomial2 = (polynomial1 / other) * other;
        result = polynomial1 - polynomial2;
        result.DeleteLeadingZeros();
        return result;
    }

    friend std::ostream& operator<<(std::ostream& out, const Poly& polynomial) {
        for (const BigInt& coefficient : polynomial.coefficients_) {
            out << coefficient << ' ';
        }
        return out;
    }

private:
    std::vector<BigInt> coefficients_;

    void DeleteLeadingZeros() {
        size_t length = 0;
        for (int64_t i = static_cast<int64_t>(coefficients_.size()) - 1; i >= 0; --i) {
            if (coefficients_[i] % p != BigInt(0)) {
                length = i + 1;
                break;
            }
        }
        coefficients_.resize(length);
    }
};

void StringToCoefficients(const std::string& str, std::vector<BigInt>& coefficients) {
    bool negative = false;
    BigInt number = 0;
    for (char i : str) {
        if (i == '-') {
            negative = true;
        } else if (i == ' ') {
            if (negative) {
                coefficients.push_back(p - number);
            } else {
                coefficients.push_back(number);
            }
            negative = false;
            number = 0;
        } else {
            number *= 10;
            number += BigInt(i - '0');
        }
    }
    if (negative) {
        coefficients.push_back(p - number);
    } else {
        coefficients.push_back(number);
    }
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


int CharToInt(char value) {
    if (value >= 48 && value <= 57) {
        return value - 48;
    }
    if (value >= 65 && value <= 90) {
        return value - 55;
    }
    if (value >= 97 && value <= 122) {
        return value - 61;
    }
    if (value == 32) {
        return 62;
    }
    if (value == 46) {
        return 63;
    }
    return 64;
}

int main() {
    std::string mess;
    std::string f;
    std::string g;
    std::string k;
    std::string message;
    std::vector<BigInt> f_polynomial;
    std::vector<BigInt> g_polynomial;
    std::vector<BigInt> k_polynomial;
    std::cin >> p;
    getline(std::cin, mess);
    getline(std::cin, f);
    getline(std::cin, g);
    getline(std::cin, k);
    getline(std::cin, message);
    StringToCoefficients(f, f_polynomial);
    StringToCoefficients(g, g_polynomial);
    StringToCoefficients(k, k_polynomial);
    Poly F(f_polynomial);
    Poly G(g_polynomial);
    Poly K(k_polynomial);
    BigInt degree = static_cast<int64_t>(f_polynomial.size()) - 1;
    std::vector<BigInt> group_message;
    BigInt number = 0;
    BigInt power = 1;
    for (char i : message) {
        number += BigInt(CharToInt(i)) * power;
        power *= 64;
    }
    while (number != 0) {
        group_message.push_back(number % p);
        number /= p;
    }
    std::vector<std::vector<BigInt>> polynomials;
    for (uint64_t i = 0; i != group_message.size(); ++i) {
        if (!(i % (f_polynomial.size() - 1))) {
            polynomials.emplace_back();
        }
        polynomials.back().push_back(group_message[i]);
    }
    for (const std::vector<BigInt>& poly : polynomials) {
        Poly message(poly);
        Poly degree_of_generator = G * G;
        degree_of_generator = degree_of_generator % F;
        std::cout << degree_of_generator << '\n';
        Poly encrypted_message = K * K;
        encrypted_message = encrypted_message * message;
        encrypted_message = encrypted_message % F;
        std::cout << encrypted_message << '\n';
    }
}



 
