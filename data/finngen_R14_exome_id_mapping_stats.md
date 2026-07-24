## Mapping Totals

| GROUP | TOTAL | ADPKD | BOTNIA | DALY | WES | PCT | NOTES |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MATCHED | 43289 | 618 | 7029 | 12223 | 23419 | 95.4% | samples with a final QRY→REF mapping in the output |
| DROPPED | 1456 | 10 | 22 | 119 | 1305 | 3.2% | found by KING but excluded from final mapping |
| NO MATCH | 625 | 0 | 110 | 47 | 468 | 1.4% | absent from ref or below KING concordance threshold |
| EXCLUDED | 5 | 0 | 0 | 0 | 5 | 0.0% | removed before KING due to sample-level QC failure |
| TOTAL | 45375 | 628 | 7161 | 12389 | 25197 | 100.0% |  |

## Mapping Breakdown

| GROUP | STATUS | TOTAL | ADPKD | BOTNIA | DALY | WES | PCT | NOTES |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MATCHED | ID_CONFIRMED | 40177 | 600 | 5936 | 11806 | 21835 | 88.5% | single candidate; KING match confirmed by matching IDs |
| MATCHED | RESOLVED_BY_ID | 160 | 6 | 31 | 102 | 21 | 0.4% | twins in ref; query ID matched one candidate |
| MATCHED | RESOLVED_BY_ALIAS | 1523 | 0 | 1053 | 0 | 470 | 3.4% | twins in ref; candidates are known aliases of each other |
| MATCHED | UNIQUE | 52 | 0 | 2 | 11 | 39 | 0.1% | single candidate; matched by genetics only |
| MATCHED | CONFLICT_KEPT | 1377 | 12 | 7 | 304 | 1054 | 3.0% | contested ref ID; kept after priority tiebreak; 1377 ref IDs contested, avg 31.4 queries/ref |
| DROPPED | CONFLICT_DROPPED | 1456 | 10 | 22 | 119 | 1305 | 3.2% | contested ref ID; lost tiebreak; REF_MAPPED = NA |
| NO MATCH | MISSING | 625 | 0 | 110 | 47 | 468 | 1.4% | no KING match found |
| EXCLUDED | HET_EXCLUDED | 5 | 0 | 0 | 0 | 5 | 0.0% | excluded by heterozygosity filter (F > 0.3) prior to KING |
