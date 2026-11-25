# tmr-software

Código e dados do artigo **"Quanto custa votar? Desempenho e cobertura de falhas da Redundância Modular Tripla em software"** (CESAR School, orientação Prof. Ronierison Maciel).

## Conteúdo

- `tmr.c`: implementação dos esquemas Base, DMR, TMR sequencial e TMR paralela, dos votadores por palavra e bit a bit, e do injetor de falhas (inversão de bit)
- `fig.py`: gera as figuras do artigo a partir dos CSVs
- `*.csv`: dados brutos usados no artigo (`tempo_mm_100reps.csv` é a repetição com 100 execuções citada nas ameaças à validade)

## Como reproduzir (Linux ou WSL2)

```bash
gcc -O2 -march=native -pthread -o tmr tmr.c
./tmr time mm 384 30        > tempo_mm.csv
./tmr time vs 16777216 30   > tempo_vs.csv
./tmr voter 1048576 30      > votador.csv
./tmr inject mm 64 1000 2027    > injecao_mm.csv
./tmr inject vs 65536 1000 2027 > injecao_vs.csv
pip install pandas scipy matplotlib
python3 fig.py   # gera fig_tempo.pdf e fig_votador.pdf
```

Ambiente usado no artigo: VM KVM, 2 vCPUs Intel Xeon 2,1 GHz, 8 GB RAM, Ubuntu 24.04, GCC 13.3.
