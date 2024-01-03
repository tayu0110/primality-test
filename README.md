# Document
Provide a method to determine whether an unsigned integer is a prime number through the `IsPrime` Trait.

In the current implementation, this crate uses the Miller-Rabin primality test.  
The Miller-Rabin primality test is known to have witnesses that can conclusively determine unsigned integers of at most 64 bits.  
In this crate, the following information is used to select the witnesses.  

[Deterministic variants of the Miller-Rabin primality test](https://miller-rabin.appspot.com/)

# LICENSE
This crate is provided under MIT License.  
Please read `LICENSE`.