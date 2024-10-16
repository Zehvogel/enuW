
```
gfortran -c Observable.f -fPIC
g++ -Wl,-soname,libObservable.so -shared -o libObservable.so Observable.o -lgfortran
nm -gDC libObservable.so
```
