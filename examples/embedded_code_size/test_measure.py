"""Check byte accounting without needing an installed MCU toolchain."""

import struct
import unittest

from measure import elf_sections, section_totals


def example_elf():
    names = b"\0.shstrtab\0.ukf_text\0.ukf_rodata\0.data\0.bss\0.debug_info\0"
    definitions = [(".shstrtab", 3, 0, names), (".ukf_text", 1, 6, b"code" * 3),
                   (".ukf_rodata", 1, 2, b"const" * 2), (".data", 1, 3, b"data"),
                   (".bss", 8, 3, b""), (".debug_info", 1, 0, b"debug" * 20)]
    payload = bytearray(52)
    entries = [bytes(40)]
    for name, kind, flags, content in definitions:
        start = len(payload)
        payload.extend(content)
        size = 1024 if kind == 8 else len(content)
        entries.append(struct.pack("<IIIIIIIIII", names.index(name.encode()), kind, flags,
                                   0, start, size, 0, 0, 1, 0))
    offset = len(payload)
    payload.extend(b"".join(entries))
    identity = b"\x7fELF\x01\x01\x01" + bytes(9)
    payload[:52] = struct.pack("<16sHHIIIIIHHHHHH", identity, 2, 40, 1, 0, 0,
                               offset, 0, 52, 0, 0, 40, len(entries), 1)
    return bytes(payload)


class SizeAccountingTests(unittest.TestCase):
    def test_code_constants_initialized_data_and_bss_are_distinct(self):
        machine, sections = elf_sections(example_elf())
        self.assertEqual(machine, 40)
        totals = section_totals(sections)
        self.assertEqual(totals["estimator_code_bytes"], 12)
        self.assertEqual(totals["estimator_rodata_bytes"], 10)
        self.assertEqual(totals["estimator_code_and_rodata_bytes"], 22)
        self.assertEqual(totals["total_flash_payload_bytes"], 26)
        self.assertEqual(totals["writable_allocated_payload_bytes"], 4)
        self.assertEqual(totals["zero_fill_bytes_including_tls"], 1024)

    def test_runtime_is_not_counted_as_estimator_code(self):
        _, sections = elf_sections(example_elf())
        sections += [{"name": ".runtime_text", "flags": 6, "type": 1, "bytes": 50},
                     {"name": ".runtime_rodata", "flags": 2, "type": 1, "bytes": 20}]
        totals = section_totals(sections)
        self.assertEqual(totals["estimator_code_and_rodata_bytes"], 22)
        self.assertEqual(totals["runtime_code_bytes"], 50)
        self.assertEqual(totals["runtime_rodata_bytes"], 20)
        self.assertEqual(totals["total_flash_payload_bytes"], 96)

    def test_rejects_wrong_class_endianness_and_truncated_sections(self):
        good = example_elf()
        bad_inputs = [b"", good[:30], good[:-1], good[:4] + b"\x02" + good[5:],
                      good[:5] + b"\x02" + good[6:]]
        for data in bad_inputs:
            with self.subTest(data=data[:10]), self.assertRaises(ValueError):
                elf_sections(data)

    def test_rejects_relocatable_objects(self):
        data = bytearray(example_elf())
        struct.pack_into("<H", data, 16, 1)
        with self.assertRaises(ValueError):
            elf_sections(data)


if __name__ == "__main__":
    unittest.main()
