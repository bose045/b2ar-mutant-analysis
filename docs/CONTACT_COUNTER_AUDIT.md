# Archived contact-frequency counter: review finding

The local GetContacts checkout inspected during the calculation-code audit reported revision
`da14deb263635e1e555519e5f208c61e83148eef`. In its `res_contacts_xl` helper:

1. A frame change flushes the preceding frame but does not append the current contact.
2. An ordinary file iterator does not yield an empty string at EOF; the final frame is not flushed.

A three-frame fixture with one identical contact in each frame returned one counted frame
(CP=1/3), while the expected result and the new counter return three frames (CP=1).
The archived frequency CLI calls this helper. This is a finding about the inspected local
code, not a claim about every upstream release or the exact engine inside an earlier container.

Reproduce the diagnostic with a trusted checkout:

```sh
python tools/audit_contact_counter.py /path/to/getcontacts/contact_calc/transformations.py
```

The new calculation path uses `analysis/contacts.py` directly and does not call that frequency
CLI. Tests verify frame transitions, EOF, duplicate atom contacts, empty frames, water endpoints
and the denominator. It also canonicalizes pairs by numeric residue indices rather than names.

## Consequences for research validation

Do not silently replace submitted CP/F values or force the new output to retain 42 edges.
To measure the actual impact, recount the same archived event files with identical windows,
classes and criteria, then compare CP, ΔCP, G, F, sign flips and retained edges. Confirm the
original execution environment first. The size/direction of manuscript-level changes has
**not** been established by the synthetic demonstration. Original research files are unchanged.

If exact historical replay is needed, retain the old frequency tables as explicitly labeled
legacy inputs. Corrected numerical results and historical reproduction are different records.
