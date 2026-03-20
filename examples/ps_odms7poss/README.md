# Example Registry

The full original runnable project has been isolated as:

- `examples/ps_odms7poss_legacy/`

That folder now contains everything from the historical layout:

- `PS-oDMS7POSS.mol`
- `PS-oDMS7POSS.pdb`
- `segment1/` to `segment4/`
- `outputs/`
- `scripts/`

Run legacy workflow from inside that example:

```bash
cd examples/ps_odms7poss_legacy/scripts
./rebuild_databases.sh
./generate_from_current_db.sh
```
