## Mapping Totals

| GROUP | TOTAL | ADPKD | BOTNIA | DALY | WES | PCT | NOTES |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MATCHED | 43200 | 624 | 7037 | 11845 | 23694 | 95.2% | samples with a final QRY→REF mapping in the output |
| DROPPED | 1558 | 5 | 16 | 507 | 1030 | 3.4% | found by KING but excluded from final mapping |
| NO MATCH | 636 | 0 | 111 | 53 | 472 | 1.4% | absent from ref or below KING concordance threshold |
| TOTAL | 45394 | 629 | 7164 | 12405 | 25196 | 100.0% |  |

## Mapping Breakdown

| GROUP | STATUS | TOTAL | ADPKD | BOTNIA | DALY | WES | PCT | NOTES |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MATCHED | ID_CONFIRMED | 28362 | 601 | 5937 | 0 | 21824 | 62.5% | single candidate; KING match confirmed by matching IDs |
| MATCHED | RESOLVED_BY_ID | 58 | 6 | 31 | 0 | 21 | 0.1% | twins in ref; query ID matched one candidate |
| MATCHED | RESOLVED_BY_ALIAS | 1529 | 0 | 1054 | 0 | 475 | 3.4% | twins in ref; candidates are known aliases of each other |
| MATCHED | UNIQUE | 11880 | 0 | 2 | 11827 | 51 | 26.2% | single candidate; matched by genetics only |
| MATCHED | CONFLICT_KEPT | 1371 | 17 | 13 | 18 | 1323 | 3.0% | contested ref ID; kept after priority tiebreak; 1371 ref IDs contested, avg 31.5 queries/ref |
| DROPPED | CONFLICT_DROPPED | 1450 | 5 | 16 | 399 | 1030 | 3.2% | contested ref ID; lost tiebreak; REF_MAPPED = NA |
| DROPPED | AMBIGUOUS_UNRESOLVED | 108 | 0 | 0 | 108 | 0 | 0.2% | multiple ref candidates; no resolution possible |
| NO MATCH | MISSING | 636 | 0 | 111 | 53 | 472 | 1.4% | no KING match found |
