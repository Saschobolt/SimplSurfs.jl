using Test
using SimplSurfs


@testset "Quad-Edge Data Structure Tests" begin

    @testset "make_edge: Pointer Consistency" begin
        mesh = PolyhedralMesh()
        e = make_edge(mesh)

        # Test rotational pointers
        @test rot(rot(rot(rot(e)))) === e
        @test rot(e) !== e
        @test rot(rot(e)) !== e

        # Test flip pointers
        @test flip(flip(e)) === e
        @test flip(e) !== e
        @test rot(rot(flip(e))) === flip(rot(rot(e))) # flip and rot commute for 180 degrees

        # Test initial 'next' pointers for an isolated edge
        @test next(e) === e
        @test next(rot(rot(e))) === rot(rot(e))
        @test next(rot(e)) === rot(rot(rot(e)))
        @test next(rot(rot(rot(e)))) === rot(e)
    end

    @testset "Data Assignment and Accessors" begin
        mesh = PolyhedralMesh()
        e = make_edge(mesh)

        v1 = Vertex(1, mesh=mesh)
        v2 = Vertex(2, mesh=mesh)
        f1 = Face(1, mesh=mesh)
        f2 = Face(2, mesh=mesh)

        # Set data using the mutator functions
        tail!(e, v1)
        head!(e, v2)
        left!(e, f1)
        right!(e, f2)

        # Test accessors
        @test tail(e) === v1
        @test head(e) === v2
        @test left(e) === f1
        @test right(e) === f2

        # Test propagation to other edges in the group
        @test tail(flip(e)) === v1 # Flipped edge shares the same tail
        @test head(flip(e)) === v2 # and head
        @test left(flip(e)) === f2 # but left/right are swapped
        @test right(flip(e)) === f1

        @test tail(rot(e)) === f2 # right(e) is tail(rot(e))
        @test head(rot(e)) === f1 # left(e) is head(rot(e))
    end

    @testset "splice!: Joining and Splitting Rings" begin
        mesh = PolyhedralMesh()
        v1 = Vertex(1, mesh=mesh)
        v2 = Vertex(2, mesh=mesh)

        a = make_edge(mesh)
        b = make_edge(mesh)
        tail!(a, v1)
        tail!(b, v2)

        # --- Test Joining ---
        # Initially, a and b are in separate 1-element rings
        @test next(a) === a
        @test next(b) === b
        @test next(flip(a)) === flip(a)
        @test next(flip(b)) === flip(b)

        # Manually unify tails to test ring joining
        # tail!(b, v1)

        # Splice with unified tails.
        SimplSurfs.splice!(a, b)

        @test next(a) === b
        @test next(b) === a
        @test next(flip(a)) === flip(b)
        @test next(flip(b)) === flip(a)


        # --- Test Splitting ---
        # Now split the 2-element ring we just created
        SimplSurfs.splice!(a, b)
        @test next(a) === a
        @test next(b) === b
        @test next(flip(a)) === flip(a)
        @test next(flip(b)) === flip(b)
    end
end



@testset "Polyhedral Mesh Tests" begin

    @testset "Polyhedral Mesh Construction" begin

        """
            _check_consistency(mesh::PolyhedralMesh)

        Check the consistency of the mesh, i.e. for each edge `e`: 
        - tail(next(e)) === tail(e) 
        - left(edge) === right(next(e))
        - tail(flip(e)) === tail(e)
        - !is_primary(e) || (!isnothing(tail(e)) && !isnothing(head(e)) && (!isnothing(left(e)) || !isnothing(right(e))))
        - !is_dual(e) || (!isnothing(left(e)) && !isnothing(right(e)) && (!isnothing(head(e)) || !isnothing(tail(e))))
        """
        function _check_consistency(mesh::PolyhedralMesh)
            for (i, e) in enumerate(edges(mesh))
                if !(tail(next(e)) === tail(e))
                    @warn "tail(next(e)) != tail(e): $(tail(next(e))) !== $(tail(e))"
                    return false
                elseif !(left(e) === right(next(e)))
                    @warn "left(e) != right(next(e)): $(left(e)) !== $(right(next(e)))"
                    return false
                elseif !(tail(flip(e)) === tail(e))
                    @warn "tail(flip(e)) != tail(e): $(tail(flip(e))) !== $(tail(e))"
                    return false
                elseif !(!is_primary(e) || (!isnothing(tail(e)) && !isnothing(head(e)) && (!isnothing(left(e)) || !isnothing(right(e)))))
                    return false
                elseif !(!is_dual(e) || (!isnothing(left(e)) && !isnothing(right(e)) && (!isnothing(head(e)) || !isnothing(tail(e)))))
                    return false
                end
            end
            return true
        end

        tetra = PolyhedralMesh([[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]])
        @test _check_consistency(tetra)
        tetra_with_hole = PolyhedralMesh([[1, 2, 3], [1, 2, 4], [1, 3, 4]])
        @test _check_consistency(tetra_with_hole)
        @test_throws ArgumentError PolyhedralMesh([[1, 2, 3], [1, 2, 4], [1, 2, 5]])
    end
end