# Parametrisierung von Pollards Rho-Methode

Dieses Repository enthält den Code für mein Jugend forscht-Projekt "Parametrisierung von Pollards Rho-Methode". Um Laufzeitmessungen durchzuführen, führe man im Projektordner

```
cargo run --release -- --measurements [M] [k_min] [k_max]
```

aus. Dabei muss `[M]` durch die Anzahl an Maschinen, `[k_min]` durch das minimale $k$ und `[k_max]` durch das maximale $k$ ersetzt werden. Es werden dann Laufzeitmessungen für alle möglichen Zuordnungen von $k$-Werten innerhalb der Schranken durchgeführt. Um die Formel auszuwerten, führe man

```
cargo run --release -- --formula [M] [k_min] [k_max]
```

aus und ersetze `[M]`, `[k_min]` und `[k_max]` entsprechend.
