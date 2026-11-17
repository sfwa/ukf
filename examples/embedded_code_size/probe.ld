/*
Link-only address layout shared by the ARM and RISC-V size probes.
These addresses are not a board memory map. No stack, heap, boot vector,
startup code, or hardware initialization is reserved or provided.

Keep probe.o separate from archive/runtime inputs to distinguish the
estimator's instantiated code and constants from required library helpers.
The driver names the compiled translation unit probe.o for every variant.
*/
ENTRY(ukf_size_initialise)

MEMORY
{
    FLASH (rx) : ORIGIN = 0x10000000, LENGTH = 64M
    RAM (rw) : ORIGIN = 0x20000000, LENGTH = 64M
}

SECTIONS
{
    .ukf_text :
    {
        *probe.o(.text .text.* .gnu.linkonce.t.*)
    } > FLASH
    .ukf_rodata :
    {
        *probe.o(.rodata .rodata.* .srodata .srodata.* .gnu.linkonce.r.*)
    } > FLASH
    .runtime_text :
    {
        *(.text .text.* .gnu.linkonce.t.*)
        *(.glue_7 .glue_7t .vfp11_veneer .v4_bx)
        *(.init .fini)
    } > FLASH
    .runtime_rodata :
    {
        *(.rodata .rodata.* .srodata .srodata.* .gnu.linkonce.r.*)
    } > FLASH
    .probe_metadata :
    {
        KEEP(*(.probe_metadata))
    } > FLASH
    .ARM.extab : { *(.ARM.extab* .gnu.linkonce.armextab.*) } > FLASH
    .ARM.exidx :
    {
        __exidx_start = .;
        *(.ARM.exidx* .gnu.linkonce.armexidx.*)
        __exidx_end = .;
    } > FLASH
    .eh_frame : { *(.eh_frame .eh_frame.*) } > FLASH
    .gcc_except_table : { *(.gcc_except_table .gcc_except_table.*) } > FLASH
    .preinit_array :
    {
        PROVIDE_HIDDEN(__preinit_array_start = .);
        KEEP(*(.preinit_array*))
        PROVIDE_HIDDEN(__preinit_array_end = .);
    } > FLASH
    .init_array :
    {
        PROVIDE_HIDDEN(__init_array_start = .);
        KEEP(*(SORT(.init_array.*))) KEEP(*(.init_array))
        PROVIDE_HIDDEN(__init_array_end = .);
    } > FLASH
    .fini_array :
    {
        PROVIDE_HIDDEN(__fini_array_start = .);
        KEEP(*(SORT(.fini_array.*))) KEEP(*(.fini_array))
        PROVIDE_HIDDEN(__fini_array_end = .);
    } > FLASH
    .data :
    {
        *(.data .data.* .gnu.linkonce.d.*)
        PROVIDE(__global_pointer$ = . + 0x800);
        *(.sdata .sdata.* .gnu.linkonce.s.*)
    } > RAM AT> FLASH
    .bss (NOLOAD) :
    {
        __bss_start__ = .;
        *(.bss .bss.* .sbss .sbss.* .gnu.linkonce.b.* .gnu.linkonce.sb.*)
        *(COMMON)
        __bss_end__ = .;
    } > RAM
    end = .;
    _end = .;
    /* Picolibc's retained allocator references these bounds. No heap is
       reserved: this ELF is linked for size inspection, never execution. */
    PROVIDE(__heap_start = end);
    PROVIDE(__heap_end = end);
    .ARM.attributes 0 : { *(.ARM.attributes) }
    .riscv.attributes 0 : { *(.riscv.attributes) }
}
