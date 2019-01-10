pt
==

[pt](https://github.com/kdbohne/pt) is an offline path tracer based on
[pbrt-v3](https://github.com/mmp/pbrt-v3), the source code provided with the
third edition of *Physically Based Rendering: From Theory to Implementation*
by Matt Phar, Wenzel Jakob, and Greg Humphreys.

This path tracer implements a subset of the functionality provided by pbrt
for the purpose of learning the intricacies of a photorealistic renderer.

Building pt
-----------

pt depends on OpenEXR. With this library installed, simply clone the repository
and run `make`.

```bash
$ git clone https://github.com/kdbohne/pt
$ cd pt
$ make
```

Running pt
----------

The build process places the `pt` executable in the `bin` directory, which can
be run with `./bin/pt` from the root of the repository. pt currently supports
`.pbrt` input files. See [pbrt-v3-scenes](https://pbrt.org/scenes-v3.html) for
example scenes.

```bash
$ ./bin/pt <input file>
```

License
-------

This source code is provided under the MIT License. Much of the code is based
on the source code of [pbrt-v3](https://github.com/mmp/pbrt-v3), which is
provided under the BSD 2-Clause "Simplified" License. The licenses of all
external dependencies are also provided. The full text of these licenses can be
found in the `LICENSE.txt` file located in this repository.
