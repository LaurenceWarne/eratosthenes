#!/usr/bin/python3

import numpy as np


def eratosthenes(limit=1000000):
    """Return a numpy array of booleans of length limit+1, where the value at
    an index corresponds to whether the index is a prime number."""
    sieve = np.ones(limit+1, dtype=bool)
    sieve[0] = sieve[1] = False
    root = int(np.ceil(limit**0.5) + 1)
    sieve[2::2] = False
    for number, boolean in enumerate(sieve[:root]):
        if boolean:
            # multiple array index assignment speeds up the
            # algorithm considerably
            sieve[number**2:limit+1:number*2] = False
    sieve[2] = True
    return sieve.nonzero()[0]


def eratosthenes_wheel(limit=1000000, wheel_size=210):
    """Return a numpy array of primes 2 to a limit inclusive, using wheel
    factorization."""
    # Get a 1 -> wheel_size array of booleans
    wheel = np.ones(wheel_size, dtype=bool)
    for p in prime_factor(wheel_size):
        wheel[p-1:wheel_size:p] = False
    q, r = divmod(limit, wheel_size)
    sieve = np.tile(wheel, q + 1)
    root = int(np.ceil(limit**0.5) + 1)
    for number, boolean in enumerate(sieve[1:root], 2):
        if (boolean):
            sieve[number**2 - 1:limit:number*2] = False
    sieve[0] = False
    sieve[np.array([prime_factor(wheel_size)]) - 1] = True
    return sieve[:-(wheel_size - r)].nonzero()[0] + 1


# TODO make a more succinct factorisation function.
def prime_factor(n):
    """Return a list of primes whose product is n."""
    prime_factors = []
    while (n & 1 == 0):
        prime_factors.append(2)
        n >>= 1
    i = 3
    while n != 1:
        if n % i == 0:
            n /= i
            prime_factors.append(i)
        else:
            i += 2
    return prime_factors
