# Security

Report a vulnerability privately through GitHub's private vulnerability reporting, at
<https://github.com/ciphernaut/dvbt-rx/security/advisories/new> (the same form as **Security**,
then **Report a vulnerability**, on the repository page; GitHub asks you to sign in). Only the
maintainer sees the report. Please do not open a public issue for it. Say what an attacker can
reach and do, with the commit and the steps to reproduce. This is a one-person project, so allow
a few days for an answer.

`dvbt-rx` reads untrusted input: IQ samples off the air, and the DVB-T signalling a transmitter
sends in them, which it turns into JSON lines. A crash, a hang or a memory error on any input is
in scope, as is anything in the stand-in headers under `shim/`. A weakness in SDRangel's
`modemdvb` library itself belongs with [its authors](https://github.com/f4exb/sdrangel).
