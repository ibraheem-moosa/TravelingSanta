import sys
import matplotlib.pyplot as plt

def get_all_primes(below):
    is_prime = [True] * below
    is_prime[0] = False
    is_prime[1] = False
    for p in range(2, below):
        if is_prime[p]:
            for q in range(2 * p, below, p):
                is_prime[q] = False
    primes = [p for p in range(below) if is_prime[p]]
    return primes

def load_cities(file_name):
    with open(file_name) as f:
        lines = [l for l in f]
        lines = [l.split(',') for l in lines]
        lines = [(int(l[0]), float(l[1]), float(l[2])) for l in lines[1:]]
        return lines

def plot_prime_cities(cities, c=None, s=1):
    primes = get_all_primes(len(cities))
    prime_cities = [cities[p] for p in primes]
    plot_cities(prime_cities, c=c, s=s)

def plot_cities(cities, c=None, s=1):
    x = [c[1] for c in cities]
    y = [c[2] for c in cities]
    plt.scatter(x, y, c=c, s=s)


if __name__ == '__main__':
    file_name = sys.argv[1]
    cities = load_cities(file_name)
    plot_cities(cities, c='b', s=0.1)
    plot_prime_cities(cities, c='r')
    plt.show()
