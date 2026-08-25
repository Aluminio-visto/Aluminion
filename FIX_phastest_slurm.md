# PHASTEST: muestras perdidas por nodos SLURM huérfanos

Commits `2c54205` + `0143ef5` · run afectado `ENTHERE_2026_JUL_07` · 2026-08-25

## Síntoma reportado

PHASTEST pasó de ~9 min/muestra a **>3 h/muestra**, justo después de
interrumpir el run (Ctrl+C) para liberar espacio en disco.

## Lo que NO era

| Hipótesis | Cómo se descartó |
|---|---|
| Disco lleno | La limpieza liberó espacio real; inodos, memoria y carga sanos |
| DB de fagos corrupta por la limpieza | `prophage_virus.db` + índices íntegros |
| Genoma difícil | Ensamblajes de muestras rápidas y lentas equivalentes |
| Caída de un nodo | Los 8 nodos aparecen en los logs |
| Contenido patológico en las piezas | Las piezas ausentes son proteína normal |
| Regresión reciente del pipeline | 2026_06_23 28/28, 2025_12_15 12/12, 2025_09_16 24/24 con `success.txt` |
| El fix de `samples` del día anterior | Muestras rápidas de esa noche ya usaban el código nuevo; los hunks del commit no tocan los módulos de plásmidos ni fagos |

## Causa raíz — dos defectos que se suman

### 1. El controlador SLURM es efímero; sus nodos no

`docker-compose.yml` **comenta** el servicio `slurmctld` independiente y da
`hostname: slurmctld` al servicio efímero que aluminion lanza con `run --rm`.
Los nodos `c1`–`c8` (`slurmd`) son contenedores **persistentes**. El cerebro del
planificador nace y muere con cada muestra mientras los nodos sobreviven, y
`slurm.conf` fija `ReturnToService=0`: un nodo que el controlador marcó caído al
morir **no vuelve nunca por sí solo**.

El propio `docker-entrypoint.sh` documenta el riesgo: deja los nodos hijo
corriendo porque es cómodo, pero *"may cause scheduling issues when phastest
container is terminated via SIGINT (Ctrl+C) before slurm jobs can finish"* —
exactamente el escenario.

**Mecanismo:** el BLAST se despacha como `sbatch --array=1-208 --wait`. Algunas
piezas del array **nunca se planifican**. Las demás acaban en ~6 min, y `--wait`
bloquea hasta que expira `phage_finder_tolerate_time=180` en `phastest.pl`.
La muestra muere sin resultados tras 3 h de espera vacía.

| Muestra | Piezas BLAST | Tiempo de BLAST | Paralelismo | Espera vacía | Piezas ausentes |
|---|---|---|---|---|---|
| SCT-HURS-43 | 207/208 | 6.3 min | 7.1× | ~173 min | 99 |
| SCT-HURS-44 | 206/208 | 6.5 min | 7.0× | ~109 min (en vivo) | 90, 208 |
| SCT-HURS-15 | 99/208 | — | 0.13× | hasta timeout | 100-208 |

Las piezas ausentes **no tienen `Job_id` en el log**: nunca se despacharon, no
es que fueran lentas. El cómputo está bien; el planificador perdió los trabajos.
La identidad de la pieza es incidental (una coincidencia inicial en la 99 hizo
pensar en determinismo; en la muestra en vuelo faltaban otras).

### 2. El centinela de resume era el directorio, que crea este propio script

`resume_done 09_phages/phastest_deep/$i` era cierto **en el instante en que la
muestra empezaba**. Una muestra interrumpida a media ejecución parecía completa
en el siguiente `--resume` y se saltaba. Sus profagos pasan a leerse como
"ninguno encontrado", **indistinguible de un negativo real**: resultado
silenciosamente incorrecto, no un fallo visible.

Observado: `SCT-HURS-13` saltada con 0 ficheros, `SCT-HURS-44` vacía y
`SCT-HURS-15/18/43` parciales (425-641 ficheros, sin `success.txt`).
Del run: **12/26 muestras** tenían `success.txt`.

> Ojo: `tRNAscan.done` es centinela de una etapa intermedia, **no** de
> finalización. Confundirlos da por completas todas las muestras.

## Fix

1. `phastest_reset_cluster()` reinicia los servicios de nodo: los re-registra con
   el siguiente controlador y limpia el estado down/drain que `ReturnToService=0`
   haría permanente. No fatal; desactivable con `$PHASTEST_NO_RESET`, lista de
   nodos en `$ALUMINION_PHASTEST_NODES`.

   > **Corregido en `0aaaa9b` — leer antes de reintroducir nada aquí.** La versión
   > original (`2c54205`) llamaba a este reset **antes de cada intento de cada
   > muestra**, es decir 26 veces en un run de 26 muestras. Eso resolvió el
   > cuelgue pero introdujo una regresión de rendimiento de 5×: además de ~13 min
   > por run en reinicios (26 × 31 s), rebotar `slurmd` contra un `slurmctld` vivo
   > deja el array de trabajos **sub-despachado**. Medido: `SCT-HURS-13`, con un
   > restart inmediatamente antes, dio ~1.2× de paralelismo efectivo sobre ocho
   > nodos (29 piezas × 15 s de CPU en 352 s de ventana) y se detuvo en seco al
   > 13% tras 17 min; `2024_AGO_20_barcode22`, sin restart previo, dio 6.8×
   > (32 piezas, 987 s de CPU en 146 s) sobre los mismos `c1`-`c8`.
   >
   > Distinción diagnóstica que importa: en `SCT-HURS-13` las piezas completadas
   > eran **1-29 contiguas** y los ocho nodos aparecían en los logs. No es un nodo
   > caído — el modo de fallo para el que se escribió este documento — sino un
   > colapso de throughput. **Ambos se ven como «PHASTEST va lento» y requieren
   > arreglos opuestos.** Por eso el reset ahora corre **una vez** al arrancar el
   > módulo (el run sí necesita ese, por los nodos que un run anterior dejó en
   > `drain`) y por muestra **sólo en el reintento**, donde ya no hay nada que
   > perder.
2. `phastest_done()` comprueba **`success.txt`**, no el directorio. Un directorio
   sin `success.txt` es residuo: se borra antes de reintentar.
3. Un reintento por muestra, condicionado a **resultado real** y no al código de
   salida — `phastest.pl` devuelve 0 tras el timeout de `phage_finder`, así que
   el exit code no distingue "sin profagos" de "no produjo nada".
4. `rm -rf JOBS/$i` **antes** del intento (`0143ef5`), no solo después: un run
   muerto entre la salida del contenedor y la limpieza deja el árbol atrás
   (`JOBS/SCT-HURS-13` como `usuario`, `JOBS/SCT-HURS-44` como `nobody:nogroup`)
   y el `cp -r` mezclaría dos intentos en una muestra.

## Validación

Sandboxes desechables: el centinela salta solo muestras con `success.txt`,
reprocesa vacías y parciales (las dos clases observadas), no salta sin
`--resume`, y limpia residuo antes de reintentar. El bucle de reintento cubre
éxito, exit-0-sin-success (timeout), exit distinto de 0, y éxito en el segundo
intento, con reinicio de nodos antes de cada uno.
`bash -n` limpio, `--help` renderiza, LF y modo 775 intactos, los 13
`scripts/*.py` parsean, **suite completa 64/64** en `aluminion_annot` (pandas 2.x).

> **Nota para quien corra la suite:** 10 tests fallan con
> `ModuleNotFoundError: No module named '_log'` si `PYTHONSAFEPATH=1` está en el
> entorno (desactiva la inserción implícita del directorio del script en
> `sys.path`, de la que depende `run_script`). Es artefacto de entorno, no
> defecto de código: `env -u PYTHONSAFEPATH python3 -m pytest tests/ -q`.

## Para relanzar

```bash
cd ~/Programs/phastest-docker && docker compose restart c1 c2 c3 c4 c5 c6 c7 c8
rm -rf ~/Programs/phastest-docker/phastest-app-docker/JOBS/SCT-HURS-{13,44}
cd /home/usuario/Seqs/Enthere
aluminion_batch --runlist runs.tsv -d ./ -b ~/Databases -t 30 -- --repo ./repository --resume
```

Se saltarán las 12 con `success.txt` y se reprocesarán 14.

## Pendiente (no implementado)

Captura de señales (`trap`) para bajar el stack limpiamente al recibir SIGINT.
El reinicio de nodos hace innecesario el arreglo manual tras cada corte, pero no
evita que el corte deje los nodos sucios en primer lugar.

**Medir si el reset inicial deprime la primera muestra.** Tras `0aaaa9b` queda
un reset por run, al arrancar el módulo. Si la muestra 1 resulta
consistentemente más lenta que las siguientes, ese reset debería volverse
condicional: resetear sólo cuando `sinfo` reporte algún nodo `down`/`drain`. Hoy
no es posible porque el `slurmctld` vive en el contenedor efímero que el pipeline
lanza por muestra, así que no hay forma de consultar el estado de SLURM desde el
host antes de arrancarlo. Resolverlo requeriría un contenedor controlador
persistente — cambio de arquitectura del stack de PHASTEST, no un parche de
`aluminion.sh`.
