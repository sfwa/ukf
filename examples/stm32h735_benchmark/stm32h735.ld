/* STM32H735IG: internal flash, AXI SRAM, and non-cacheable DTCM only. */
ENTRY(Reset_Handler)
MEMORY
{
    FLASH (rx) : ORIGIN = 0x08000000, LENGTH = 1024K
    AXIRAM (xrw) : ORIGIN = 0x24000000, LENGTH = 320K
    DTCMRAM (xrw) : ORIGIN = 0x20000000, LENGTH = 128K
}
_estack = ORIGIN(DTCMRAM) + LENGTH(DTCMRAM);
__stack_limit = _estack - 64K;
SECTIONS
{
    .isr_vector : { . = ALIGN(4); KEEP(*(.isr_vector)) } >FLASH
    .text :
    {
        . = ALIGN(4);
        *(.text*) *(.glue_7) *(.glue_7t)
        KEEP(*(.init)) KEEP(*(.fini))
        . = ALIGN(4);
    } >FLASH
    .rodata : { . = ALIGN(4); *(.rodata*) . = ALIGN(4); } >FLASH
    .ARM.extab : { *(.ARM.extab* .gnu.linkonce.armextab.*) } >FLASH
    .ARM.exidx :
    {
        __exidx_start = .; *(.ARM.exidx*) __exidx_end = .;
    } >FLASH
    .preinit_array :
    {
        PROVIDE_HIDDEN(__preinit_array_start = .);
        KEEP(*(.preinit_array*))
        PROVIDE_HIDDEN(__preinit_array_end = .);
    } >FLASH
    .init_array :
    {
        PROVIDE_HIDDEN(__init_array_start = .);
        KEEP(*(SORT(.init_array.*))) KEEP(*(.init_array*))
        PROVIDE_HIDDEN(__init_array_end = .);
    } >FLASH
    .fini_array :
    {
        PROVIDE_HIDDEN(__fini_array_start = .);
        KEEP(*(SORT(.fini_array.*))) KEEP(*(.fini_array*))
        PROVIDE_HIDDEN(__fini_array_end = .);
    } >FLASH
    _sidata = LOADADDR(.data);
    .data :
    {
        . = ALIGN(4); _sdata = .; *(.data*) . = ALIGN(4); _edata = .;
    } >AXIRAM AT>FLASH
    .bss (NOLOAD) :
    {
        . = ALIGN(4); _sbss = .; __bss_start__ = .;
        *(.bss*) *(COMMON)
        . = ALIGN(4); _ebss = .; __bss_end__ = .;
    } >AXIRAM
    .heap (NOLOAD) :
    {
        . = ALIGN(8); __heap_start = .; end = .; _end = .;
        . += 16K; __heap_end = .;
    } >AXIRAM
    /* Debugger reads bypass the CPU cache. DTCM keeps the capture coherent. */
    .capture (NOLOAD) :
    {
        . = ALIGN(8); KEEP(*(.capture)) . = ALIGN(8);
    } >DTCMRAM
    ASSERT(ADDR(.capture) + SIZEOF(.capture) <= __stack_limit, "Capture overlaps benchmark stack")
    .stack __stack_limit (NOLOAD) : { . += 64K; } >DTCMRAM
    .ARM.attributes 0 : { *(.ARM.attributes) }
}
